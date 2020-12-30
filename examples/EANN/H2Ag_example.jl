
#push!(LOAD_PATH, pwd())
using NonadiabaticMolecularDynamics.IO
using NonadiabaticMolecularDynamics.Systems
using NonadiabaticMolecularDynamics.Models.EANN_H2Ag
using Unitful
using UnitfulAtomic

model_path = "/Users/wojciechstark/Desktop/ML-model-repository-master/H2_on_Ag111/h2ag111pes"
input_f = model_path * "/Ag.xyz"
atoms, positions = read_system(input_f)
model = EannH2AgModel(model_path, atoms)
p = Systems.System(atoms, model)
positions_ang = copy(ustrip(auconvert.(u"Å", positions)))

println("Initialize...")
EANN_H2Ag.initialize_H2Ag_pes(model_path)
println("Positions:")
println(positions_ang)
println("Energy:")
en = model.get_V0(positions_ang)
println(en)
println("Forces:")
f = model.get_D0(positions_ang)
println(f)


