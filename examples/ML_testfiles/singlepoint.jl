using NonadiabaticMolecularDynamics.IO
using NonadiabaticMolecularDynamics.Models.ML
using PyCall
ase = pyimport("ase")
model_path = "examples/ML_testfiles"
atoms, phasespace = 
read_system("examples/ML_testfiles/H2.xyz")
model = SchNetPackModel(model_path, atoms)
# model 1 : H2 in gas pahse
R = ase.io.read(joinpath(model_path,"H2.xyz")).get_positions()
R= reshape(R,length(R))
energy = model.get_V0(R)
forces = model.get_D0(R)
print(energy)
print(forces)


# model 1 : H2 in gas pahse
atoms2, phasespace2 = read_system("examples/ML_testfiles/H2_Cu(100).xyz")
model2 = SchNetPackModel(model_path, atoms2)
R2 = ase.io.read(joinpath(model_path,"H2_Cu(100).xyz")).get_positions()
R2= reshape(R2,length(R2))
energy2 = model2.get_V0(R2)
print(energy2)