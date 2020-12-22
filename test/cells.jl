using Test
using NonadiabaticMolecularDynamics

# 1D cells
@test PeriodicCell(hcat(1)) isa PeriodicCell
@test PeriodicCell(hcat(1.0)) isa PeriodicCell

# 3D cells
x = [1.0, 0.0, 0.0]
y = [0.0, 1.0, 0.0]
z = [0.0, 0.0, 1.0]
@test PeriodicCell([x y z]) isa PeriodicCell

a = PeriodicCell(hcat(1))
@test a.periodicity == [true, true, true]
set_periodicity!(a, [true, true, false])
@test a.periodicity == [true, true, false]

# InfiniteCell
@test InfiniteCell{Float64}() isa InfiniteCell