using NonadiabaticModels
using NonadiabaticMolecularDynamics
using Random, Distributions

using PyCall
using Unitful, UnitfulAtomic

"""
    Purpose: Integrate Belal's implementation of NO/Au(111) so we
    can run our respective IESH implementations against each other and debug 
    them together.
"""

Random.seed!(13)
# Number of electronic states
M = 4
nstates = M



# # # Test James's model (problem here is that Ase seems to handle the positions in some non-intutive way)
# # build = pyimport("ase.build")
# # io = pyimport("ase.io")
# # # #slab = build.fcc111("Au", size=(11,12,4), vacuum=10.0, a=4.175)
# # slab = build.fcc111("Au", size=(2,2,4), vacuum=10.0, a=4.175)
# # no = build.molecule("NO")
# # build.add_adsorbate(slab, no, 1.5, "ontop")
# # #slab = io.read("/home/sjanke/Documents/Uni/Warwick/Documents/Projects/Anderson_Holstein_H/IESH_Tully/My_IESH/data/surface_Au111.in")
# # temperature = austrip(300u"K")
# # atoms, r, cell = NonadiabaticDynamicsBase.convert_from_ase_atoms(slab)
# # println(r)
# # v = zero(r)
# # println(v)


# 1b) Initialize atomic parameters, i.e., moving atoms, which is NO in this case
#atoms = Atoms{Float64}([:N, :O])

# Set up geometries for calculation
Å = 1.0e-10 # 1 Angström in meter
timestep = 0.1e-15 #fs
bohr = 0.529177210903
#println(auconvert(timestep*u"s"))

filename = "surface_Au111_4.dat"
#filename = "/home/sjanke/Documents/Uni/Warwick/Documents/Projects/Anderson_Holstein_H/IESH_Tully/My_IESH/data/surface_Au111.dat"
f = open(filename)
data = readlines(f)
close(f)
n_au = parse(Int64, data[1])
a0 = parse(Float64, data[2])*sqrt(2)*Å
cell_mat = zeros(3)
cell_mat[1] = parse(Float64, data[3])*Å
cell_mat[2] = parse(Float64, data[4])*Å
cell_mat[3] = parse(Float64, data[2])*sqrt(6)*100*Å
au_atoms = zeros(n_au,3)
no_pos = zeros(2,3)

# N and O 1.15 A appart
no_pos[1,3] = 10.118647440105047#*Å
no_pos[2,3] = 10.881352559894953#*Å

#no_pos[1,3] = 1.5*Å
#no_pos[2,3] = 2.65*Å
no_pos[1,3] = 1.5*bohr
no_pos[2,3] = 2.65*bohr

for i in 5:length(data)
    au_atoms[i-4,:] = parse.(Float64,split(strip(data[i])))#*Å
end

atoms=Atoms(vcat([:N, :O], fill(:Au, Int(length(au_atoms)/3))))

# Get new model.
p = zeros(n_au+2,3)
p[1:2,:] = no_pos
#p[3:n_au+2,:] = au_atoms*Å
p[3:n_au+2,:] = au_atoms*bohr

model = Tully_NOAu111(M = M+1, n_au=n_au, au_pos=au_atoms, no_pos = no_pos,
                      x_pos = p, cell_mat = cell_mat, a_lat = a0, W=7.0)



#Still needs: positions handed into potential and derivative
# p are the positions, I believe
pto = potential(model, p)
dto = derivative(model, p)

#println(size(dto))
#println(size(dto[1]))
#println(dto[1])
#println(pto)
#@test Dynamics.IESH{Float64}(42,40) isa Dynamics.IESH
# # #Initialize the simulation problem; Simulation is defined in src/simulation_constructors.jl 

sim = Simulation{IESH_Tully}(atoms, model, n_electrons=Int(M/2))

# From comparison with James's implemementation, r = (x1,x2,...,xn; y1,y2,...,yn; z1,z2,..., zn)
# v should have the same format
#v = fill(5/sim.atoms.masses[1], sim.DoFs, length(sim.atoms))
r1 = p'
v = zero(r1)
v[3,1:2] .= 0.1
v = ustrip.(map( x -> auconvert(x*u"Å/fs"), v))

# Check forces
#for I in eachindex(v)
#    println(dto[I][1,1:2])
#end
#r = fill(-5.0, sim.DoFs, length(sim.atoms))

z = DynamicsVariables(sim,v, r1)

@time solution = run_trajectory(z, (0.0, 30auconvert(timestep*u"s")), sim, dt=auconvert(timestep*u"s"), adaptive=false; 
                                            output=(:position))
# println("Finished")

plot(solution, :position)
