using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Models
using PyCall

#os=pyimport("os")
#os.environ['KMP_DUPLICATE_LIB_OK']='True'
#ase = pyimport("ase") # is this needed?
model_path = "examples/ML_testfiles"
frictionmodel_path = "examples/ML_testfiles/friction/"

# model 1 : H2 in gas phase
cell, atoms, positions = read_system("examples/ML_testfiles/H2.xyz")
model = Models.SchNetPackModels.SchNetPackModel(model_path, cell, atoms)

V = [0.0]
potential!(model, V, positions)
display(V)

D = zero(positions)
derivative!(model, D, positions)
display(D)

#test for friction model

# ! Attention ! 
# Install friction_tensor model (https://github.com/mgastegger/friction_tensor) 
# You need permission from Michael Gastegger.
# you also need to install oyaml in pip

#for two atoms manually. could this be changed somehow to the size 3*natoms times 3*natoms? but 
F = Any
#TODO dimension of F should be 6x6
cell_friction, atoms_friction, positions_friction = read_system("examples/ML_testfiles/H2.xyz")
frictionmodel = Models.SchNetPackModels.FrictionSchNetPackModel(frictionmodel_path, cell_friction, atoms_friction)
#frictionmodel = Models.SchNetPackModels.SchNetPackModel(frictionmodel_path, cell_friction, atoms_friction)

friction!(frictionmodel, F, positions_friction)
display(F)

#D = zero(positions)
#derivative!(model, D, positions)
#display(D)
