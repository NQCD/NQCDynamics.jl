using Test
using NQCBase

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
@test a.tmp_vector1 == [0, 0, 0]
@test a.tmp_vector2 == [0, 0, 0]
@test a.tmp_bools == [false, false, false]

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
@test InfiniteCell() isa InfiniteCell

@testset "evaluate_periodic_distance" begin
    cell = PeriodicCell([1 0 0; 0 1 0; 0 0 1])
    r1 = [0.1, 0.1, 0.1]
    r2 = [0.9, 0.9, 0.9]
    @test sqrt(3*0.2^2) ≈ evaluate_periodic_distance(cell, r1, r2)
    @test sqrt(3*0.2^2) ≈ evaluate_periodic_distance(cell, r2, r1)
    r1 = [0.4, 0.4, 0.4]
    r2 = [0.6, 0.6, 0.6]
    @test sqrt(3*0.2^2) ≈ evaluate_periodic_distance(cell, r1, r2)
    @test sqrt(3*0.2^2) ≈ evaluate_periodic_distance(cell, r2, r1)
end

@testset "check_atoms_in_cell" begin
    cell = PeriodicCell([1 0 0; 0 1 0; 0 0 1])
    R = rand(3, 10)
    @test check_atoms_in_cell(cell, R) # All atoms inside cell
    R[3, 5] += 1
    @test !check_atoms_in_cell(cell, R) # Some atoms not inside cell
    R[3, 5] -= 2
    @test !check_atoms_in_cell(cell, R) # Some atoms not inside cell
end

@testset "Periodic translations" begin
    cell = PeriodicCell([1 0 0; 0 1 0; 0 0 1])
    R = rand(3, 10)
    translation_xy = NQCBase.PeriodicReplica(cell, [1,1,0])
    R_xy = translation_xy(R)
    # Does the translation work as expected?
    for idx in axes(R_xy, 2)
        @test R_xy[:, idx] == R[:, idx] + [1,1,0]
    end
    # Try building a Supercell which contains a single translation
    supercell_xy = Supercell(cell, 1,1,0) 
    R_sup = supercell_xy(R)
    # The supercell only contains the requested replicas, the initial positions are not prepended.
    @test R_sup == R_xy
    # Replicas in all directions include the [0,0,0] translation exactly once
    supercell_3x3 = Supercell(cell, -1:1, -1:1, 0)
    R_3x3 = supercell_3x3(R)
    @test size(R_3x3) == (3, 9 * size(R, 2))
    # Iterators.product runs x fastest, so [0,0,0] is the 5th replica
    @test R_3x3[:, 41:50] == R
    @test count(j -> R_3x3[:, j] == R[:, 1], axes(R_3x3, 2)) == 1
end

@testset "Supercell of a Structure" begin
    cell = PeriodicCell([1 0 0; 0 1 0; 0 0 1])
    R = rand(3, 2)
    structure = NQCBase.Structure(Atoms([:H, :C]), R, cell)
    supercell = Supercell(cell, -1:1, 0, 0)
    replicated = supercell(structure)
    @test replicated isa NQCBase.Structure
    @test replicated.positions == supercell(R)
    @test replicated.atoms.types == [:H, :C, :H, :C, :H, :C]
    @test size(replicated.positions, 2) == length(replicated.atoms.types)
    @test replicated.cell === structure.cell
end