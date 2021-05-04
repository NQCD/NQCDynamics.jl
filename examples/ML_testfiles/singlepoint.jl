using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Models
using PyCall
#os=pyimport("os")
#os.environ['KMP_DUPLICATE_LIB_OK']='True'
#ase = pyimport("ase") # is this needed?
model_path = "examples/ML_testfiles"

# model 1 : H2 in gas phase
cell, atoms, positions = read_system("examples/ML_testfiles/H2.xyz")
model = SchNetPackModel(model_path, cell, atoms)

V = [0.0]
potential!(model, V, positions)
display(V)

D = zero(positions)
derivative!(model, D, positions)
display(D)

# model 1 : H2 in gas phase
#cell, atoms, positions = read_system("examples/ML_testfiles/H2_Cu(100).xyz")
#model = SchNetPackModel(model_path, cell, atoms)##

#V = [0.0]
#potential!(model, V, positions)
#display(V)

#D = zero(positions)
#derivative!(model, D, positions)
#display(D)
