using NQCModels
using NQCDynamics
# using TullyNOAu111
using FortranNOAu111
using Random
using LinearAlgebra
using Unitful, UnitfulAtomic
using ProfileView
using DelimitedFiles
using Libdl: DL_LOAD_PATH
push!(DL_LOAD_PATH, "/home/sjanke/Documents/Uni/Warwick/Documents/Projects/Anderson_Holstein_H/fork/FortranNOAu111/lib")

Random.seed!(13)
# Number of electronic states
#M = 40
M = 40
nstates = M


Å = 1.0e-10 # 1 Angström in meter
timestep = 0.5u"fs"
bohr = 0.529177210903

#filename = "./test/Dynamics/surface_Au111.dat"
filename = "/home/sjanke/Documents/Uni/Warwick/Documents/Projects/Anderson_Holstein_H/IESH_Tully/My_IESH/data/surface_Au111.dat"
#filename = "surface_Au111.dat"
#filename = "surface_Au111_l1a.dat"
f = open(filename)
data = readlines(f)
close(f)
n_au = parse(Int64, data[1])
a0 = parse(Float64, data[2])*sqrt(2)*Å
cell_mat = zeros(3)
cell_mat[1] = parse(Float64, data[3])*Å
cell_mat[2] = parse(Float64, data[4])*Å
cell_mat[3] = parse(Float64, data[2])*sqrt(6)*100*Å
au_atoms = zeros(n_au,3)
no_pos = zeros(2,3)

# N and O 1.15 A appart
#no_pos[1,3] = ustrip(auconvert(1.5*u"Å"))
#no_pos[2,3] = ustrip(auconvert(2.65*u"Å"))
n = 7
no_pos[1,1] = ustrip(auconvert(15.269558583375413*u"Å"))
no_pos[1,2] = ustrip(auconvert(3.1520123873416158*u"Å"))
no_pos[1,3] = ustrip(auconvert((11.569538179345007)*u"Å"))
no_pos[2,1] = ustrip(auconvert(15.325962733588399*u"Å"))
no_pos[2,2] = ustrip(auconvert(4.0472432519126691*u"Å"))
no_pos[2,3] = ustrip(auconvert((10.848669353650543)*u"Å"))

for i in 5:length(data)
    au_atoms[i-4,:] = parse.(Float64,split(strip(data[i])))#*Å
end

atoms=Atoms(vcat([:O, :N], fill(:Au, Int(length(au_atoms)/3))))

p = zeros(n_au+2,3)
p[1:2,:] = no_pos
p[3:n_au+2,:] = map( x -> ustrip(auconvert(x*u"Å")), au_atoms)

#model = FortranNOAu111Model("../FortranNOAu111.jl/lib/tullynoau111", permutedims(p); Ms=M)
#model = FortranNOAu111Model("/home/sjanke/Documents/Uni/Warwick/Documents/Projects/Anderson_Holstein_H/fork/FortranNOAu111/lib/tullynoau111", permutedims(p); Ms=M, freeze_layers=4)
model = FortranNOAu111Model(permutedims(p); Ms=M, freeze_layers=4)

r1 = permutedims(p)

sim = Simulation{DiabaticIESH}(atoms, model, n_electrons=Int(M/2))

v = zero(r1)
v[3,1:2] .= -2000
v = ustrip.(map( x -> auconvert(x*u"m/s"), v))

z = DynamicsVariables(sim,v, r1)

# Get Energies
#DynamicsUtils.classical_potential_energy(sim,z)
# Get forces
#println(map( x -> ustrip(auconvert(u"J/m", x)), NQCDynamics.DynamicsOutputs.force(z, 0.0, sim)[:,1:3])')

@time solution = run_trajectory(z, (0.0, 2*timestep), sim, dt=timestep; 
                                            output=(:position, :hamiltonian, :force))


write_extxyz("traj.xyz", atoms, solution.position, PeriodicCell(diagm(cell_mat)))


##print results to file
open("trajectory.txt", "w") do fi
    write(fi, "#step, Nx, Ox, Ny, Oy, Nz, Oz\n")
    writedlm(fi, [(ustrip.(auconvert.(u"m", p))[:,1:3])' for p in solution.position], '\t')
    #writedlm(fi, [(ustrip.(auconvert.(u"m", p))[:,1:2])' for p in solution.position], '\t')
end

open("energy.txt", "w") do fi
    write(fi, "#step, energy (eV)\n")
    writedlm(fi, [(ustrip.(auconvert.(u"J", p)))' for p in solution.hamiltonian], '\t')
end

open("forces.txt", "w") do fi
    write(fi, "#step, Nx, Ox, Ny, Oy, Nz, Oz\n")
    writedlm(fi, [(ustrip.(auconvert.(u"J/m", p))[:,1:3])' for p in solution.force], '\t')
end

###############################################################################
#
# My Script version
#
###############################################################################
# using NQCModels
# using NQCDynamics
# using TullyNOAu111
# using Random, Distributions
# using LinearAlgebra
# using FiniteDiff
# using Test

# using PyCall
# using Unitful, UnitfulAtomic
# using ProfileView
# using DelimitedFiles

# """
#     Purpose: Integrate Belal's implementation of NO/Au(111) so we
#     can run our respective IESH implementations against each other and debug 
#     them together.
# """

# # function finite_difference_gradient(model::NQCModels.DiabaticModels.DiabaticModel, R)
# #     f(x, i, j) = potential(model, x)[i,j]
# #     grad = [Hermitian(zeros(NQCModels.nstates(model), NQCModels.nstates(model))) for i in CartesianIndices(R)]
# #     for k in eachindex(R)
# #         for i=1:NQCModels.nstates(model)
# #             for j=1:NQCModels.nstates(model)
# #                 grad[k].data[i,j] = FiniteDiff.finite_difference_gradient(x->f(x,i,j), R)[k]
# #             end
# #         end
# #     end
# #     grad
# # end

# function finite_difference_gradient(model::NQCModels.DiabaticModels.DiabaticModel, R)
#     f(x, i, j) = potential(model, x)[i,j]
#     grad = [Hermitian(zeros(NQCModels.nstates(model), NQCModels.nstates(model))) for i in CartesianIndices(R)]
#     for i=1:NQCModels.nstates(model)
#         for j=1:NQCModels.nstates(model)
#             gradient = FiniteDiff.finite_difference_gradient(x->f(x,i,j), R)
#             for k in eachindex(R)
#                 grad[k].data[i,j] = gradient[k]
#             end
#         end
#     end
#     grad
# end

# Random.seed!(13)
# # Number of electronic states
# #M = 40
# M = 2
# nstates = M



# # # # Test James's model (problem here is that Ase seems to handle the positions in some non-intutive way)
# # # build = pyimport("ase.build")
# # # io = pyimport("ase.io")
# # # # #slab = build.fcc111("Au", size=(11,12,4), vacuum=10.0, a=4.175)
# # # slab = build.fcc111("Au", size=(2,2,4), vacuum=10.0, a=4.175)
# # # no = build.molecule("NO")
# # # build.add_adsorbate(slab, no, 1.5, "ontop")
# # # #slab = io.read("/home/sjanke/Documents/Uni/Warwick/Documents/Projects/Anderson_Holstein_H/IESH_Tully/My_IESH/data/surface_Au111.in")
# # # temperature = austrip(300u"K")
# # # atoms, r, cell = NonadiabaticDynamicsBase.convert_from_ase_atoms(slab)
# # # println(r)
# # # v = zero(r)
# # # println(v)


# # 1b) Initialize atomic parameters, i.e., moving atoms, which is NO in this case
# #atoms = Atoms{Float64}([:N, :O])

# # Set up geometries for calculation
# Å = 1.0e-10 # 1 Angström in meter
# #timestep = 0.1e-15 #fs
# timestep = 0.5u"fs"
# bohr = 0.529177210903
# #println(auconvert(timestep*u"s"))

# filename = "/home/sjanke/Documents/Uni/Warwick/Documents/Projects/Anderson_Holstein_H/IESH_Tully/My_IESH/data/surface_Au111_4.dat"
# #filename = "/home/sjanke/Documents/Uni/Warwick/Documents/Projects/Anderson_Holstein_H/IESH_Tully/My_IESH/data/surface_Au111.dat"
# #filename = "surface_Au111.dat"
# filename = "surface_Au111_l1.dat"
# f = open(filename)
# data = readlines(f)
# close(f)
# n_au = parse(Int64, data[1])
# a0 = parse(Float64, data[2])*sqrt(2)*Å
# cell_mat = zeros(3)
# cell_mat[1] = parse(Float64, data[3])*Å
# cell_mat[2] = parse(Float64, data[4])*Å
# cell_mat[3] = parse(Float64, data[2])*sqrt(6)*100*Å
# au_atoms = zeros(n_au,3)
# no_pos = zeros(2,3)

# # N and O 1.15 A appart
# #no_pos[1,3] = ustrip(auconvert(1.5*u"Å"))
# #no_pos[2,3] = ustrip(auconvert(2.65*u"Å"))
# n = 7
# no_pos[1,1] = ustrip(auconvert(15.269558583375413*u"Å"))
# no_pos[1,2] = ustrip(auconvert(3.1520123873416158*u"Å"))
# no_pos[1,3] = ustrip(auconvert((11.569538179345007-n)*u"Å"))
# no_pos[2,1] = ustrip(auconvert(15.325962733588399*u"Å"))
# no_pos[2,2] = ustrip(auconvert(4.0472432519126691*u"Å"))
# no_pos[2,3] = ustrip(auconvert((10.848669353650543-n)*u"Å"))

# for i in 5:length(data)
#     au_atoms[i-4,:] = parse.(Float64,split(strip(data[i])))#*Å
# end

# atoms=Atoms(vcat([:O, :N], fill(:Au, Int(length(au_atoms)/3))))

# # Get new model.
# p = zeros(n_au+2,3)
# p[1:2,:] = no_pos
# #p[3:n_au+2,:] = au_atoms*Å
# #map( x -> auconvert(u"J", x), V)

# p[3:n_au+2,:] = map( x -> ustrip(auconvert(x*u"Å")), au_atoms)



# model = Tully_NOAu111(M = M+1, n_au=n_au, au_pos=au_atoms, no_pos = no_pos,
#                       x_pos = p, cell_mat = cell_mat, a_lat = a0, W=7.0)



# #Still needs: positions handed into potential and derivative
# # p are the positions, I believe
# r1 = p'
# pto = potential(model, r1)
# dto = derivative(model, r1)
# println(map( x -> ustrip(auconvert(u"J", x)), pto))

# #println(size(dto))
# #println(size(dto[1]))
# #println(dto[1])
# #println(pto)
# #@test Dynamics.IESH{Float64}(42,40) isa Dynamics.IESH
# # # #Initialize the simulation problem; Simulation is defined in src/simulation_constructors.jl 

# #sim = Simulation{IESH_Tully}(atoms, model, n_electrons=Int(M/2))
# sim = Simulation{DiabaticIESH}(atoms, model, n_electrons=Int(M/2))

# # From comparison with James's implemementation, r = (x1,x2,...,xn; y1,y2,...,yn; z1,z2,..., zn)
# # v should have the same format
# #v = fill(5/sim.atoms.masses[1], sim.DoFs, length(sim.atoms))
# v = zero(r1)
# #v[3,1:2] .= -0.5
# #x, N
# # v[1,1] = -46.024809328916717
# # v[2,1] = -730.49287493313341
# # v[3,1] = -1948.0774185178323
# # v[1,2] = 40.293357512232333
# # v[2,2] = 639.52487797328479
# # v[3,2] = -3051.2599916426143
# v[3,1:2] .= -2000
# #v = ustrip.(map( x -> auconvert(x*u"Å/fs"), v))
# v = ustrip.(map( x -> auconvert(x*u"m/s"), v))


# z = DynamicsVariables(sim,v, r1)

# @time solution = run_trajectory(z, (0.0, 10*timestep), sim, dt=timestep, adaptive=false; 
#                                             output=(:position))


# @time solution = run_trajectory(z, (0.0, 600*timestep), sim, dt=timestep, adaptive=false; 
#                                             output=(:position))



# ##print results to file
# open("trajectory.txt", "w") do fi
#     write(fi, "#step, Nx, Ox, Ny, Oy, Nz, Oz\n")
#     writedlm(fi, [(ustrip.(auconvert.(u"m", p))[:,1:2])' for p in solution.position], '\t')
# end


# #println([(ustrip.(auconvert.(u"m", p))[:,1:2])' for p in solution.position])

# # Test finite differences and compare to derivatives.
# #function test_model(model::NQCModels.Model, R; rtol=1e-5)
# # rtol=1e-5
# # D = derivative(model, r1)
# # println(D)
# # finite_diff = finite_difference_gradient(model, r1)
# # println(finite_diff)
# # isapprox(finite_diff, D, rtol=rtol)