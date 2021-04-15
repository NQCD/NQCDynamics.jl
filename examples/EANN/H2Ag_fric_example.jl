using NonadiabaticMolecularDynamics
using PyCall

model_path = "~/programs/ML-model/repository/H2_on_Ag111"
input = model_path * "h2ag111pes/Ag.xyz"
cell, atoms, positions = read_system(input)

println("Initialize...")
model = Models.EANN_H₂Ag(model_path, atoms)

println("Friction:")
friction = zeros(3*length(atoms),3*length(atoms))
Models.friction!(model, friction, positions)

println("friction: / 1/ps")
display(friction[1:10,1:10])
