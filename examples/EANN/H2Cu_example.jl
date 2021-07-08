using NonadiabaticMolecularDynamics
using PyCall

model_path = "/Users/wojciechstark/Desktop/H2_on_Cu/1_h2cu_pes"
input_f = model_path * "/H2Cu_example.xyz"
cell, atoms, positions = read_system(input_f)

println("Initialize...")
model = EANN_H₂Cu(model_path, atoms)

println("Positions:")
println(positions) 

println("Energy:")
V = [0.0]
potential!(model, V, positions)
println(V)

println("Forces:")
D = zero(positions)
derivative!(model, D, positions)
println(D)

println("Deallocate...")
deallocate_H2Cu_pes(model_path)
