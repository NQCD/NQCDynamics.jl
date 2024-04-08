using Test
using NQCDynamics
using Unitful, UnitfulAtomic
using PyCall
using NQCDynamics, NQCModels
using JLD2

@testset "Desorption Analysis" begin
    ase_io=pyimport("ase.io")
    desorption_trajectory_ase=ase_io.read("artifacts/desorption_test.xyz", index=":")
    desorption_dynamicsvariables=jldopen("artifacts/desorption_dynamics.jld2")["trajectory"]
    # Make a dummmy simulation out of the first structure (we don't need a working potential to test. )
    atoms, initial_positions, cell=NQCDynamics.NQCBase.convert_from_ase_atoms(desorption_trajectory_ase[1])
    simulation=Simulation(atoms, AdiabaticASEModel(desorption_trajectory_ase[1]), cell=cell)
    diatomic_indices=[55,56]
    @test @time "DesorptionFrame - No views" get_desorption_frame_noviews(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=2.4u"Å") == 2675
    @test @time "DesorptionFrame - With views" get_desorption_frame_views(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=2.4u"Å") == 2675
    @test @time "DesorptionAngle - No views" get_desorption_angle(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=2.4u"Å", evaluation_function=get_desorption_frame_noviews) ≈ 28.05088202518
    @test @time "DesorptionAngle - With views" get_desorption_angle(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=2.4u"Å", evaluation_function=get_desorption_frame_views) ≈ 28.05088202518
end

