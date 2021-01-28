using NonadiabaticMolecularDynamics
using PyCall

model_path = "/Users/wojciechstark/Desktop/ML-model-repository-master/H2_on_Ag111/h2ag111pes"
input_f = model_path * "/Ag.xyz"
cell, atoms, positions = read_system(input_f)

println("Initialize...")
model = Models.EANN_H₂Ag(model_path, atoms)

println("Positions:")
println(positions)

println("Energy:")
V = [0.0]
Models.potential!(model, V, positions)
println(V)

println("Forces:")
D = zero(positions)
Models.derivative!(model, D, positions)
println(D)
