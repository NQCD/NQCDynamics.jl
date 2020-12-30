using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO

model_path = "/Users/wojciechstark/Desktop/ML-model-repository-master/H2_on_Ag111/h2ag111pes"
input_f = model_path * "/Ag.xyz"
cell, atoms, positions = read_system(input_f)

println("Initialize...")
model = Models.EannH2AgModel(model_path, atoms)

println("Positions:")
println(positions)

println("Energy:")
V = [0.0]
model.potential!(V, positions)
println(V)

println("Forces:")
D = zero(positions)
model.derivative!(D, positions)
println(D)
