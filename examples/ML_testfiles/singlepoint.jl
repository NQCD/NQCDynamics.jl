using NonadiabaticMolecularDynamics.IO
using NonadiabaticMolecularDynamics.Models.ML
using PyCall
ase = pyimport("ase")
model_path = "examples/ML_testfiles"

# model 1 : H2 in gas pahse
cell, atoms, positions = read_system("examples/ML_testfiles/H2.xyz")
model = SchNetPackModel(model_path, cell, atoms)

V = [0.0]
model.potential!(V, positions)
display(V)

D = zero(positions)
model.derivative!(D, positions)
display(D)

# model 1 : H2 in gas pahse
cell2, atoms2, phasespace2 = read_system("examples/ML_testfiles/H2_Cu(100).xyz")
model2 = SchNetPackModel(model_path, cell2, atoms2)

V = [0.0]
model.potential!(V, positions)
display(V)

D = zero(positions)
model.derivative!(D, positions)
display(D)
