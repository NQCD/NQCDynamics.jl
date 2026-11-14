
using Test
using NQCBase

atoms = Atoms([:H, :C, :O, :N])
cell = PeriodicCell(rand(3, 3) .* 10)
R = rand(3, 4) .* 10
structure = NQCBase.Structure(atoms, R, cell)

@testset "to/from_extxyz_dict (Atoms, Positions, Cell)" begin
    dict = NQCBase.to_extxyz_dict(atoms, R, cell) # Test conversion with atoms, R, cell and Structure methods. 
    @test dict["cell"] ≈ au_to_ang.(permutedims(cell.vectors, (2,1)))
    converted_structure = NQCBase.from_extxyz_dict(dict)
    @test converted_structure.cell.vectors ≈ cell.vectors
    @test converted_structure.cell.inverse ≈ cell.inverse
end
@testset "to/from_extxyz_dict (Structure)" begin
    dict = NQCBase.to_extxyz_dict(structure) # Test conversion with atoms, R, cell and Structure methods. 
    @test dict["cell"] ≈ au_to_ang.(permutedims(structure.cell.vectors, (2,1)))
    converted_structure = NQCBase.from_extxyz_dict(dict)
    @test converted_structure.cell.vectors ≈ structure.cell.vectors
    @test converted_structure.cell.inverse ≈ structure.cell.inverse
end

@testset "Single frame" begin
    file_buffer = "output.xyz"
    write_extxyz(file_buffer, atoms, R, cell) # Save a 1-structure file
    structure = read_extxyz(file_buffer) |> first # Load the 1-structure file
    @test structure.atoms == atoms
        @test all(isapprox.(structure.cell.vectors, cell.vectors, atol = 1e-8))
        @test all(isapprox.(structure.cell.inverse, cell.inverse; atol = 1e-8))
        @test structure.cell.periodicity ≈ cell.periodicity
    @test all(isapprox.(structure.positions, R, atol = 1e-8))
end

@testset "Multiple frames" begin
    fb = "output.xyz"
    atoms = Atoms([:H, :C, :O, :N])
    cell = PeriodicCell(rand(3, 3) .* 10)
    R = [rand(3, 4) .* 10 for _=1:100]
    write_extxyz(fb, atoms, R, cell)
    structures = read_extxyz(fb)
    for (i,structure) in enumerate(structures)
        @test structure.atoms == atoms
        @test all(isapprox.(structure.cell.vectors, cell.vectors, atol = 1e-8))
        @test all(isapprox.(structure.cell.inverse, cell.inverse; atol = 1e-8))
        @test structure.cell.periodicity ≈ cell.periodicity
        @test all(isapprox.(structure.positions, R[i], atol = 1e-8))
    end
end

@testset "Non-periodic structures" begin
    @testset "Single non-periodic frame" begin
        file_buffer = "non_periodic_single.xyz"
        atoms = Atoms([:H, :C, :O, :N])
        cell = InfiniteCell()
        R = rand(3, 4) .* 10
        structure_in = NQCBase.Structure(atoms, R, cell)

        # Write and read back
        write_extxyz(file_buffer, structure_in)
        structure_out = read_extxyz(file_buffer) |> first

        # Verify atoms and positions
        @test structure_out.atoms == atoms
        @test all(isapprox.(structure_out.positions, R, atol = 1e-8))

        # Verify cell is InfiniteCell
        @test structure_out.cell isa InfiniteCell

        rm(file_buffer)
    end

    @testset "Multiple non-periodic frames" begin
        file_buffer = "non_periodic_multiple.xyz"
        atoms = Atoms([:H, :C, :O, :N])
        cell = InfiniteCell()
        R = [rand(3, 4) .* 10 for _=1:50]
        structures_in = [NQCBase.Structure(atoms, R[i], cell) for i=1:50]

        # Write and read back
        write_extxyz(file_buffer, structures_in)
        structures_out = read_extxyz(file_buffer)

        # Verify each structure
        for (i, structure) in enumerate(structures_out)
            @test structure.atoms == atoms
            @test all(isapprox.(structure.positions, R[i], atol = 1e-8))
            @test structure.cell isa InfiniteCell
        end

        rm(file_buffer)
    end

    @testset "Write non-periodic with Atoms, R, Cell" begin
        file_buffer = "non_periodic_arc.xyz"
        atoms = Atoms([:H, :C, :O])
        cell = InfiniteCell()
        R = rand(3, 3) .* 5

        # Write using the atoms, R, cell interface
        write_extxyz(file_buffer, atoms, R, cell)
        structure = read_extxyz(file_buffer) |> first

        @test structure.atoms == atoms
        @test all(isapprox.(structure.positions, R, atol = 1e-8))
        @test structure.cell isa InfiniteCell

        rm(file_buffer)
    end

    @testset "Multiple non-periodic with Atoms, R, Cell" begin
        file_buffer = "non_periodic_arc_multi.xyz"
        atoms = Atoms([:H, :C, :O])
        cell = InfiniteCell()
        R = [rand(3, 3) .* 5 for _=1:25]

        # Write using the atoms, R, cell interface
        write_extxyz(file_buffer, atoms, R, cell)
        structures = read_extxyz(file_buffer)

        for (i, structure) in enumerate(structures)
            @test structure.atoms == atoms
            @test all(isapprox.(structure.positions, R[i], atol = 1e-8))
            @test structure.cell isa InfiniteCell
        end

        rm(file_buffer)
    end
end

rm("output.xyz")
