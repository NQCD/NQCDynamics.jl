using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO

input_f = model_path * "/H2Cu_example.xyz"
cell, atoms, positions = read_system(input_f)

println("Initialize...")
model = Models.EANN_H₂Cu(model_path, atoms)

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
EANN_H2Cu.deallocate_H2Cu_pes(model_path)
