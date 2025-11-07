using Test
using NQCDynamics
using LinearAlgebra

@testset "Surface site classification" begin
   test_cell = PeriodicCell(I(3) .* 3) # Unit length cell
   dummy_sites = Dict(
    :type_1 => [
        [0.0,0.0],
    ],
    :type_2 => [
        [1.0,1.0],
    ],
   )
   testing_positions = [rand(3,10) for _ in 1:4]
   # Indices 9 and 10 are the atom to test
   testing_positions[1][:, 9] .= [0,0,0]
   testing_positions[1][:, 10] .= [0,0,0]
   testing_positions[2][:, 9] .= [1,1,0]
   testing_positions[2][:, 10] .= [1,1,0]
   testing_positions[3][:, 9] .= [0.5,0,0]
   testing_positions[3][:, 10] .= [0.5,0,0]
   testing_positions[4][:, 9] .= [0.25,0,0]
   testing_positions[4][:, 10] .= [0.25,0,0]
   # Case 1: 0,0 should snap to site type 1
   result_case1 = Analysis.HighSymmetrySites.positions_to_category(
    testing_positions[1][:, 9],
    dummy_sites,
    PeriodicCell(I(2))
   )
   @test result_case1 == :type_1
   # Case 2: 1,1 should snap to site type 2
   result_case2 = Analysis.HighSymmetrySites.positions_to_category(
    testing_positions[2][:, 9],
    dummy_sites,
    PeriodicCell(I(2))
   )
   @test result_case2 == :type_2
   # Case 3: 0.5,0 should map to :other if snap_to_site isn't particularly large. 
   result_case3 = Analysis.HighSymmetrySites.positions_to_category(
    testing_positions[3][:, 9],
    dummy_sites,
    PeriodicCell(I(2))
   )
   @test result_case3== :other
   # Case 4: 0.25,0 should map to site type 1 if snap_to_site is set to 0.5. 
   result_case4 = Analysis.HighSymmetrySites.positions_to_category(
    testing_positions[2][:, 9],
    dummy_sites,
    PeriodicCell(I(2)),
    snap_to_site = 0.5,
   )
   @test result_case4 == :type_1
   # Case 5: Analysing the full trajectory should give the correct results. 
   result_case5 = Analysis.HighSymmetrySites.classify_every_frame(
    [DynamicsVariables(v = rand(3,10), r = i) for i in testing_positions],
    test_cell,
    Analysis.HighSymmetrySites.SlabStructure(
        [9,10],
        dummy_sites,
        [3.0,3.0,1.0],
    ),
    snap_to_site = 0.45,
   )
   for i in 1:2
    @test result_case5[1][i] == :site_1
    @test result_case5[2][i] == :site_2
    @test result_case5[3][i] == :other
    @test result_case5[4][i] == :site_1
   end
end