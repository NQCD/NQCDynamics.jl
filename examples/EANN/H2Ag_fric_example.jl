using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO

model_path = "/Users/wojciechstark/Desktop/ML-model-repository-master/H2_on_Ag111/h2ag111pes"
input_f = model_path * "/Ag_t.xyz"
input1 = model_path * "/Ag.xyz"
cell, atoms, positions = read_system(input1)

println("Initialize...")
model = Models.EANN_H₂Ag(model_path, atoms)

println("Friction:")
open(input_f) do file
    for ln in eachline(file)
        m = split(ln)
        if size(m)[1] > 5 && m[1]=="H"
            friction = zeros(6,6)
            coor = zeros(6,1)
            for i in 1:6
                coor[i,1]=parse(Float64, m[i+1])
            end
            coor_friction = copy(transpose(coor))
            Models.friction!(model, friction, coor_friction)
            println(friction)
        end
    end
end

