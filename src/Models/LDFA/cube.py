#! /usr/bin/env python
import numpy as np
import sys
import time

def arr2D(list2D):
    '''arr2D converts a 2D list into a two-dimensional np.array. This is faster than converting a 2D list with the function np.array,
    because the latter involves a lot of checks on the dimensions of the elements of the list.'''
    h = len(list2D)
    w = len(list2D[0])
    array = np.empty((h,w))
    for i in range(h):
        for j in range(w):
            array[i,j] = list2D[i][j]
    return array

def arr1D(list1D):
    length = len(list1D)
    array = np.empty(length)
    for i in range(length):
        array[i] = list1D[i]
    return array

class cube:
    '''The cube class represents an object used for storing and manipulating data read from cube files. '''
    def __init__(self):
        self.bohr2ang = 0.52917721092
        self.n_atoms = None
        self.origin = None
        self.x_len = None
        self.y_len = None
        self.z_len = None
        self.n_points = None
        self.x_vec = None
        self.y_vec = None
        self.z_vec = None
        self.xyz_array = None
        self.jacobi = None
        self.filename = None
        self.atoms_array = None
        self.density = None
        self.cube_file_set = False

    def __call__(self,pos_x, pos_y,pos_z, silent=False):
        """
        This command takes three real numbers that correspond to a
        Cartesian x, y, and z position and returns the value of the 
        Voxel that is closest to this position.
        """

        vec = [pos_x,pos_y,pos_z]
        if self.cube_file_set:
            cube_dimensions = self.xyz_array*[self.x_len,self.y_len,self.z_len]
            cube_dim_inv = np.linalg.inv(cube_dimensions)
            frac_coord = np.dot(cube_dim_inv.T,vec-self.origin)
            #frac_coord = np.dot(vec-self.origin,cube_dim_inv)
            if any([i>1.00 for i in frac_coord]) or any([i<0.0 for i in frac_coord]):
                if silent:
                    pass
                else:
                    print('point lies outside of cube cell = density 0.0')
                return 0.0
            int_coord = np.array(frac_coord * [self.x_len,self.y_len,self.z_len],dtype=np.int)
            copy_density = self.density.reshape([self.x_len,self.y_len,self.z_len])
            return copy_density[int_coord[0],int_coord[1],int_coord[2]]
        else:
            print('First set a cube file with self.read')

    def read(self, filename, castep2cube_format=False):
        ''' The function reads a cube file starting from the third line (the first two lines are usually reserved for comments and contain no data).
        cube.read returns a one-dimensional np.array in which each element resembles a value of the density data in the cube file.
        The function also sets the variables: n_atoms, origin, x_len, y_len, z_len, n_points, x_vec, y_vec, z_vec, xyz_array, density and jacobi
        of the corresponding cube object.'''

        if filename:
            self.filename=str(filename)
        cubefile = open(self.filename, "r")
        print(("Reading cube file: ", cubefile.name))
        next(cubefile)
        next(cubefile)
        cube_content_list = [line for line in cubefile]
        cube_content_joined = "".join(cube_content_list)
        cubefile.close()
        cube_content_split = cube_content_joined.split()
        cube_content = list(map(float, cube_content_split))
        self.n_atoms = int(cube_content[0])
        self.origin = np.array(cube_content[1:4])*self.bohr2ang
        self.x_len = int(cube_content[4])
        self.y_len = int(cube_content[8])
        self.z_len = int(cube_content[12])
        self.n_points = self.x_len*self.y_len*self.z_len
        self.x_vec = np.array(cube_content[5:8])
        self.y_vec = np.array(cube_content[9:12])
        self.z_vec = np.array(cube_content[13:16])
        self.xyz_array = np.array([self.x_vec,self.y_vec,self.z_vec])*self.bohr2ang
        self.jacobi = np.linalg.det(self.xyz_array)
        self.atoms_array = arr1D(cube_content[16:16+5*self.n_atoms])
        self.atoms_array.shape = (self.n_atoms, 5)
        self.cube_file_set = True
        if castep2cube_format:
            self.density = arr1D(cube_content[16+5*self.n_atoms:len(cube_content)]) / self.jacobi
        else:
            self.density = arr1D(cube_content[16+5*self.n_atoms:len(cube_content)])

    def write(self, filename_out):
        '''writes out a read-in cube file with correct formatting of the data'''
        fd=open('./'+str(filename_out), 'w+')
        fd.write('Density:\n')
        fd.write('\n')
        fd.write(str('%5i %12.6f %12.6f %12.6f\n') %(self.n_atoms, self.origin[0]/self.bohr2ang, self.origin[1]/self.bohr2ang, 
            self.origin[2]/self.bohr2ang))
        fd.write(str('%5i %12.6f %12.6f %12.6f\n') %(self.x_len, self.x_vec[0], self.x_vec[1], self.x_vec[2]))
        fd.write(str('%5i %12.6f %12.6f %12.6f\n') %(self.y_len, self.y_vec[0], self.y_vec[1], self.y_vec[2]))
        fd.write(str('%5i %12.6f %12.6f %12.6f\n') %(self.z_len, self.z_vec[0], self.z_vec[1], self.z_vec[2]))
        for i in range(len(self.atoms_array)):
                fd.write(str('%5i') %(self.atoms_array[i][0]))
                fd.write(str('%12.6f %12.6f %12.6f') %(self.atoms_array[i][1], self.atoms_array[i][2], \
                                                       self.atoms_array[i][3]))
                fd.write(str('%12.6f') %(self.atoms_array[i][4]))
                fd.write('\n')
        for i in range(len(self.density)):
            fd.write(str('%12.6f') %(self.density[i]))
            if (i%6==5):
                fd.write('\n')
        fd.close()

    def scale(A,a=1.0):
        new_density = A.density*a
        A.density=new_density

    def add_and_scale(self,A,a=1.0):
        new_density = self.density + a*A.density
        self.density = new_density

    def __add__(A,B):
        '''Adds the density of a cube object B to A cube object A.
        This is not necessarily a commutative operation,
        as A possibly contains other atoms than B.'''
        new_density = A.density + B.density
        A.density = new_density

    def __sub__(A,B):
        new_density = A.density - B.density
        A.density = new_density

    def crop(self):
        self.density.shape = (self.x_len,self.y_len,self.z_len)
        self.density = self.density[0:self.x_len-1,0:self.y_len-1,0:self.z_len-1]
        self.density = self.density.flatten()
        self.x_len = self.x_len-1
        self.y_len = self.y_len-1
        self.z_len = self.z_len-1
        print((type(self.density)))
        print((self.density.shape))

    def integrate_xy(self, filename_out):
        copy_density = self.density
        copy_density.shape = (self.x_len,self.y_len,self.z_len)
        copy_density = copy_density.sum(axis=(0,1))*np.linalg.norm(np.cross(self.x_vec,self.y_vec))*self.bohr2ang**2
        print(('shape of original density', self.density.shape))
        print(('shape of new density', copy_density.shape))
        fd=open('./'+str(filename_out), 'w+')
        for i in range(len(copy_density)):
            fd.write(str('%12.6f %12.6f\n') %( i*np.linalg.norm(self.z_vec)*self.bohr2ang,copy_density[i]))
        fd.close()

    def dipole_of_cores(self):
        array = np.array(self.atoms_array)
        core_charges = {'1':1,'6':4,'7':5,'29':11,'47':11}
        dipole = 0
        for core in array:
            dipole += core_charges[str(int(core[0]))]*core[2:]*self.bohr2ang
        return dipole

    def multipole(self, order=int(0)):
        '''multipole returns the multipole moment as np.array of the requested
        order (an integer) related to the specified cube objects density in np.array format
        and the volume of a voxel in the cube file (jacobi).'''

        print('Calculating multipole moment')
        if order>0:
            ind_array = [[integer/self.z_len/self.y_len,integer/self.z_len%self.y_len,integer%self.z_len] \
                         for integer in range(self.n_points)]
            index_array = arr2D(ind_array)
            r_vec = np.dot(index_array,self.xyz_array)
        if order == 0:
            result = np.sum(self.density)*self.jacobi
            print(('multipole moment of order 0:\n', result))
            return result
        elif order>0:
            outer = np.outer
            def multiple_outer(a):
                b = a
                for i in range(order-1):
                    b = outer(b,a)
                b.shape=tuple([3]*order)
                return b
            mult_array=list(map(multiple_outer, r_vec))
        else:
            raise ValueError('Value for the input parameter "order" is not defined.')
        mult_array = np.array(mult_array)
        result = np.tensordot(self.density, mult_array,axes=1)*self.jacobi
        result.shape = tuple([3]*order)
        print(('multipole moment of order %i:\n' %order, result))
        return result

    def dipole_map_xy(self):
        '''calculates an xy map of the dipole moments in z direction'''

        print('Calculating dipole moment map')
        e0 = 0.00552635

        copy_density = self.density
        copy_density.shape = (self.x_len,self.y_len,self.z_len)
        r_vec = [float(i*self.xyz_array[2,2]) for i in range(self.z_len)]
        result = np.zeros([self.x_len,self.y_len])
        for x in range(self.x_len):
            for y in range(self.y_len):

                den = copy_density[x,y,:]
                result[x,y] = np.dot(den,r_vec)#* self.xyz_array[2,2]

        workfunction = result / e0 / np.linalg.norm(np.cross(self.x_len*self.xyz_array[0],\
                                                             self.y_len*self.xyz_array[1]))
        return result, workfunction

