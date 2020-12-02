push!(LOAD_PATH,pwd())
push!(LOAD_PATH,joinpath(pwd(),"src","Dynamics"))

#this is just for testing and copied from the examples of MDPrePostProcessing
using NonadiabaticMolecularDynamics.Systems
using NonadiabaticMolecularDynamics.Models
using NonadiabaticMolecularDynamics.Models.ML
using NonadiabaticMolecularDynamics.Models.ML_descriptor
#include(joinpath(pwd(),"src","Models","ML","ML.jl"))
using Unitful
using UnitfulAtomic
using Revise

using PeriodicTable


#test run
model_path=joinpath(string(pwd()),"examples","ML_testfiles")

x = [1.0, 0.0, 0.0]
y = [0.0, 1.0, 0.0]
z = [0.0, 0.0, 1.0]

a = PeriodicCell([x y z])
set_periodicity!(a, [false, true, false])
p = AtomicParameters(PeriodicCell(zeros(3, 3)), [:H, :H])
system = System(p, Models.Analytic.Free())
R=[5.08124009 0.94499003 8.42028042; 4.72585008 0.70419003 7.57972028]
#init ML
property_dict = Dict("energy" => true, "forces" => true)
model,model_args,force_mask = initialize_MLmodel(model_path,p,a)
#init schnet_descriptor
#environment_provider=env_init(model_args)

@time begin
schnet_input=julia2schnet(p,R,model_args)
schnet_outputs=model(schnet_input)
#60ms

end
@time begin
t=1000
macro run_prediction(n,exp)
    quote
        for i = 1:$(esc(n))
            $(esc(exp))
        end
    end
end

#@run_prediction(t,(model(schnet_input)))
end
print(schnet_outputs)
print(schnet_input)