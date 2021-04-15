using NonadiabaticMolecularDynamics
using PyCall
using ASE
using Unitful
using Plots
import JuLIP

build = pyimport("ase.build")
ase = pyimport("ase")

slab = build.fcc111("Pd", size=(3,2,3), orthogonal=true, vacuum=10.0)
build.add_adsorbate(slab, "H", 1.5, "ontop")
build.add_adsorbate(slab, "H", 1.5, "bridge")

cell, atoms, positions = extract_parameters_and_positions(slab)
model = Models.JuLIPModel(atoms, cell, JuLIP.Potentials.EAM("PdH_hijazi.eam.alloy"))

Δ = Dict([(:H, 5.0)])
sim = Simulation(atoms, model; temperature=10000u"K", cell=cell)

"""
Make sure the hydrogens do not enter the surface
"""
function constrain_above_surface(Rₚ)
    (au_to_ang(Rₚ[3,end]) > 15) && (au_to_ang(Rₚ[3,end-1]) > 15)
end

output = InitialConditions.run_monte_carlo_sampling(sim, positions, Δ, 100;
                                fix=collect(1:18), extra_function=constrain_above_surface)

@show output.acceptance
write_trajectory("sampling.xyz", cell, atoms, output.R)
plot(output.energy)