using NonadiabaticMolecularDynamics.Atoms
using NonadiabaticMolecularDynamics.Systems
using NonadiabaticMolecularDynamics.Models
using NonadiabaticMolecularDynamics.Models.ML
using NonadiabaticMolecularDynamics.Models.ML_descriptor
using Unitful
using UnitfulAtomic
using PeriodicTable
using BenchmarkTools

model_path=joinpath(string(pwd()),"examples","ML_testfiles")

x = [1.0, 0.0, 0.0]
y = [0.0, 1.0, 0.0]
z = [0.0, 0.0, 1.0]

a = PeriodicCell([x y z])
set_periodicity!(a, [false, true, false])
p = AtomicParameters(a, [:H, :H])
R = [5.08124009 0.94499003 8.42028042; 4.72585008 0.70419003 7.57972028]
R = [5.08124009, 0.94499003, 8.42028042, 4.72585008, 0.70419003, 7.57972028]

property_dict = Dict("energy" => true, "forces" => true)
model, model_args, force_mask = initialize_MLmodel(model_path, p)

@btime julia2schnet($p, $R, $model_args) # 751 μs, 450 allocations
inputs = Dict{String, PyObject}()
@btime update_schnet_input!($inputs, $p, $R, $model_args) # 727 μs, 391 allocations

@btime $model($inputs) # 5.78 ms, 47 allocations