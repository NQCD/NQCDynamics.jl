using Test
using NonadiabaticMolecularDynamics

# 1D cells
@test PeriodicCell(hcat(1)) isa PeriodicCell
@test PeriodicCell(hcat(1.0)) isa PeriodicCell

# 3D cells
x = [1.0, 0.0, 0.0]
y = [0.0, 1.0, 0.0]
z = [0.0, 0.0, 1.0]
a = PeriodicCell([x y z])
@test a isa PeriodicCell
@test a.inverse == a.vectors

@testset "apply_cell_boundaries!" begin
    cell = PeriodicCell([1 0 0; 0 1 0; 0 0 1])
    R = rand(3, 4)
    A = copy(R)
    apply_cell_boundaries!(cell, A)
    @test R == A # Check unchanged when inside cell
    A += 2rand(3, 4) # Move atoms out of cell
    apply_cell_boundaries!(cell, A)
    @test all(0 .<= A .<= 1) # Check they're all back in
    A -= 2rand(3, 4) # Move atoms out of cell
    apply_cell_boundaries!(cell, A)
    @test all(0 .<= A .<= 1) # Check they're all back in
end

a = PeriodicCell(hcat(1))
@test a.periodicity == [true, true, true]
set_periodicity!(a, [true, true, false])
@test a.periodicity == [true, true, false]

# InfiniteCell
@test InfiniteCell{Float64}() isa InfiniteCell