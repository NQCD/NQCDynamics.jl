#test run
model_path=joinpath(string(pwd()),"examples","ML_testfiles")
#this is just for testing and copied from the examples of MDPrePostProcessing
using NonadiabaticMolecularDynamics.Systems
using NonadiabaticMolecularDynamics.Models.ML
using Unitful
using PyCall
using UnitfulAtomic
using PeriodicTable
ase = pyimport("ase")
spk = pyimport("schnetpack")
x = [0.0, 0.0, 0.0]
y = [0.0, 0.0, 0.0]
z = [0.0, 0.0, 0.0]
filename="H2.xyz"
a = PeriodicCell([x y z])
set_periodicity!(a, [false, false, false])
p = AtomicParameters(PeriodicCell(zeros(3, 3)), [:H, :H])
system = System(p, Models.Analytic.Free())
R = rand(6)
print(R)
#read file and get atoms object
device="cpu"
atoms = ase.io.read(joinpath(model_path,filename),index=-1,format="xyz")
environment_provider = spk.utils.script_utils.settings.get_environment_provider(model_args,device=device)
py"""
def pyf(atoms):
    print(atoms.pbc.astype("uint8"))
"""
py"pyf"(atoms)
#init ML
property_dict = Dict("energy" => true, "forces" => true)
model, model_args, force_mask = initialize_MLmodel(model_path,p,a)
calculator = spk.interfaces.SpkCalculator(model = model, device=device, energy = spk.Properties.energy, forces=spk.Properties.forces,energy_units="eV", forces_units="eV/A",environment_provider=environment_provider)
@time begin
atoms = ase.io.read(joinpath(model_path,filename),index=-1,format="xyz")
print(atoms.get_positions())
end
@time begin
#-4.187468528747559Float32[-1.383376 -0.93732774 -3.2719274; 1.383376 0.93732774 3.2719274]  2.449919 seconds (4.54 M allocations: 206.210 MiB, 3.89% gc time)
atoms.set_calculator(calculator)
energy = atoms.get_total_energy()
forces=atoms.get_forces()
end
print(energy,forces)
#predict
#force mask will be either an array or false
#first array of tuple = energy, second array of tuple = forces

#running dynamics
#how long
t=10
macro run_prediction(n,exp)
    quote
        for i = 1:$(esc(n))
            $(esc(exp))
        end
    end
end
#@run_prediction(t,println(get_properties(calculator,atoms)))
#0.130670 seconds (33.89 k allocations: 1.696 MiB)