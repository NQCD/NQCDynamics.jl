using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO

model_path = "/Users/wojciechstark/Desktop/ML-model-repository-master/NO_on_Au111"
input_f = model_path * "/NOAu_example.xyz"
cell, atoms, positions = read_system(input_f)

println("Initialize...")
model = Models.EANN_NOAu(model_path, atoms)

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

println("Deallocate...")
Models.deallocate_NOAu_pes(model_path)
