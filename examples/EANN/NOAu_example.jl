using NonadiabaticMolecularDynamics
using PyCall

model_path = "/Users/wojciechstark/Desktop/ML-model-repository-master/NO_on_Au111"
input_f = model_path * "/NOAu_example.xyz"
cell, atoms, positions = read_system(input_f)

println("Initialize...")
model = Models.EANN_NOAu(model_path, atoms)

println("Energy:")
V = [0.0]
Models.potential!(model, V, positions)
println(V)

println("Forces:")
D = zero(copy(transpose(positions)))
println(D)
Models.derivative!(model, D, positions)
println(D)

println("Friction")
F = zeros(6,6)
Models.friction!(model, F, positions)
print(F)

println("Deallocate...")
Models.deallocate_NOAu_pes(model_path)
