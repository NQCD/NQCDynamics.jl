using Test
using NQCDynamics
using Unitful, UnitfulAtomic
using PythonCall
using NQCModels
using JLD2

@testset "Desorption Analysis" begin
    ase_io=pyimport("ase.io")
    desorption_trajectory_ase=ase_io.read("artifacts/desorption_test.xyz", index=":")
    desorption_dynamicsvariables=jldopen("artifacts/desorption_dynamics.jld2")["trajectory"]
    # Make a dummmy simulation out of the first structure (we don't need a working potential to test. )
    atoms, initial_positions, cell=NQCDynamics.convert_from_ase_atoms(desorption_trajectory_ase[1])
    simulation=Simulation(atoms, AdiabaticASEModel(desorption_trajectory_ase[1]), cell=cell)
    diatomic_indices=[55,56]
    @test Analysis.Diatomic.get_desorption_frame(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=austrip(2.4u"Å")) == 2696
    @test Analysis.Diatomic.get_desorption_angle(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=austrip(2.4u"Å")) ≈ 25.445786265092522
end
