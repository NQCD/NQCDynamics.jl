
#push!(LOAD_PATH, pwd())
using NonadiabaticMolecularDynamics.IO
using NonadiabaticMolecularDynamics.Systems
using NonadiabaticMolecularDynamics.Models.EANN_H2Cu
using Unitful
using UnitfulAtomic

model_path = "/Users/wojciechstark/Desktop/H2_on_Cu/1_h2cu_pes/"
input_f = model_path * "/H2Cu_example.xyz"
atoms, positions = read_system(input_f)
model = EANN_H2Cu.EannH2CuModel(model_path, atoms)
p = Systems.System(atoms, model)

println("Initialize...")
EANN_H2Cu.initialize_H2Cu_pes(model_path)
println("Positions:")
println(positions) 
println("Energy:")
en = model.get_V0(positions)
println(en)
println("Forces:")
f = model.get_D0(positions)
println(f)
println("Deallocate...")
EANN_H2Cu.deallocate_H2Cu_pes(model_path)
