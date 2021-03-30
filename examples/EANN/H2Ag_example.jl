using NonadiabaticMolecularDynamics
using PyCall

pes_path = "/home/chem/msuhmw/programs/ML-model-repository/H2_on_Ag111/h2ag111pes/h2pes.so"
friction_path = "/home/chem/msuhmw/programs/ML-model-repository/H2_on_Ag111/h2agtensor/tensor_neu/h2agtensor.so"
input_f = "Ag.xyz"
cell, atoms, positions = read_system(input_f)

println("Initialize...")
model = Models.EANN_H₂Ag(pes_path, friction_path, atoms)

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
