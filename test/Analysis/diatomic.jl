using Test
using PythonCall
using NQCDynamics
using Unitful, UnitfulAtomic
using NQCModels
using NQCBase
using JLD2

@testset "Desorption Analysis" begin
    ase_io=pyimport("ase.io")
    desorption_trajectory_ase=ase_io.read("artifacts/desorption_test.xyz", index=":")
    desorption_dynamicsvariables=jldopen("artifacts/desorption_dynamics.jld2")["trajectory"]
    # Make a dummmy simulation out of the first structure (we don't need a working potential to test. )
    atoms, initial_positions, cell=convert_from_ase_atoms(desorption_trajectory_ase[1])
    simulation=Simulation(atoms, AdiabaticASEModel(desorption_trajectory_ase[1]), cell=cell)
    diatomic_indices=[55,56]
    @test Analysis.Diatomic.get_desorption_frame(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=2.4u"Å") == 2675
    @test Analysis.Diatomic.get_desorption_angle(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=2.4u"Å") ≈ 28.05088202518
end