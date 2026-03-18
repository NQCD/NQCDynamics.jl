using Test
using NQCDynamics
using Unitful, UnitfulAtomic
using PythonCall
using NQCModels
using NQCDInterfASE
using LinearAlgebra
using JLD2

@testset "Desorption Analysis" begin
    ase_io=pyimport("ase.io")
    desorption_trajectory_ase=ase_io.read("artifacts/desorption_test.xyz", index=":")
    desorption_dynamicsvariables=jldopen("artifacts/desorption_dynamics.jld2")["trajectory"]
    # Make a dummy simulation out of the first structure (we don't need a working potential to test. )
    structure=convert_from_ase_atoms(desorption_trajectory_ase[1])
    simulation=Simulation(structure.atoms, ClassicalASEModel(desorption_trajectory_ase[1]), cell=structure.cell)
    diatomic_indices=[55,56]
    @test Analysis.Diatomic.get_desorption_frame(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=austrip(2.4u"Å")) == 2675
    @test Analysis.Diatomic.get_desorption_angle(desorption_dynamicsvariables, diatomic_indices, simulation; surface_distance_threshold=austrip(2.4u"Å")) ≈ 28.05088202518
end

@testset "Internal coordinate transformation" begin
    structure_import = NQCDynamics.read_extxyz("artifacts/desorption_test.xyz")
    # Take only the first structure
    first_structure = NQCDynamics.NQCBase.Structure(structure_import[1], structure_import[2][1], structure_import[3])
    # Generate a random (symmetric) friction tensor
    eft = Symmetric(rand(6,6)) |> Matrix
    diatomic_indices = [55,56]
    simulation=Simulation(first_structure.atoms, Free(), cell=first_structure.cell)
    # Ensure internal coordinate transformation is trace-conserving
    eft_internal = NQCDynamics.Analysis.Diatomic.transform_to_internal_coordinates(
        eft,
        first_structure.positions,
        diatomic_indices...,
        simulation
    )
    @test isapprox(eft |> diag |> sum, eft_internal |> diag |> sum)
    # Ensure doing the back-transform on the forward transform gives back the same matrix. 
    eft_cartesian_again = NQCDynamics.Analysis.Diatomic.transform_from_internal_coordinates(
        eft_internal,
        first_structure.positions,
        diatomic_indices...,
        simulation
    )
    matrix_comparison = isapprox.(eft, eft_cartesian_again)
    for idx in eachindex(matrix_comparison)
        @test matrix_comparison[idx]
    end
end