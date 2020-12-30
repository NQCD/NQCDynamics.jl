using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO

model_path = "/Users/wojciechstark/Desktop/H2_on_Cu/1_h2cu_pes"
input_f = model_path * "/H2Cu_example.xyz"
cell, atoms, positions = read_system(input_f)

println("Initialize...")
model = Models.EANN_H₂Cu(model_path, atoms)

println("Positions:")
println(positions) 

println("Energy:")
V = zeros(1, 1)
model.potential!(V, positions)
println(V)

println("Forces:")
D = zero(positions)
model.derivative!(D, positions)
println(D)

println("Deallocate...")
Models.deallocate_H2Cu_pes(model_path)
