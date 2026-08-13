using RingPolymerArrays
using Test
using LinearAlgebra: I, SymTridiagonal, Diagonal
using Statistics: mean

nbeads = 2

A = [
    1 4
    2 5
    3 6
]
ndofs, natoms = size(A)

@testset "IndexMapping" begin
    classical = [2, 4, 8]
    quantum = RingPolymerArrays.find_quantum_indices(10, classical)
    map = RingPolymerArrays.IndexMapping(quantum, classical)
    @test map.classical[2] == 1
    @test map.classical[4] == 2
    @test map.classical[8] == 3
    @test map.quantum[3] == 2
    @test map.quantum[5] == 3
end

@testset "constructor" begin
    classical = [2, 4]
    rp = RingPolymerArray{Float64}(undef, (3, 5, 3); classical)
    @test collect(RingPolymerArrays.classicalindices(rp)) == [2, 4]
    @test collect(RingPolymerArrays.quantumindices(rp)) == [1, 3, 5]
end

@testset "getindex" begin
    rp = RingPolymerArray(A; nbeads)
    for I in CartesianIndices(rp)
        @test rp[I] == A[I[1], I[2]]
    end
end

@testset "setindex!" begin
    rp = RingPolymerArray(A; nbeads)
    rp[3,1,2] = 0
    @test rp.quantum_atoms[2][3,1] == 0
    rp = RingPolymerArray(A; nbeads, classical=[1])
    rp[3,1,1] = 0
    @test rp.classical_atoms[3,1] == 0
end

@testset "similar" begin
    rp = RingPolymerArray(A; nbeads)
    rp_similar = similar(rp)
    @test rp.classical_atoms !== rp_similar.classical_atoms
    @test rp.quantum_atoms !== rp_similar.classical_atoms
    @test rp.index_map.classical == rp_similar.index_map.classical
    @test rp.index_map.quantum == rp_similar.index_map.quantum
    @test rp.dims == rp_similar.dims
end

@testset "eachbead" begin
    A = rand(ndofs, natoms, nbeads)
    rp = RingPolymerArray(A)
    for (i,bead) in enumerate(eachbead(rp))
        @test bead == A[:,:,i]
    end
end

@testset "eachdof" begin
    A = rand(ndofs, natoms, nbeads)
    rp = RingPolymerArray(A)
    for (I, dof) in zip(CartesianIndices(A[:,:,1]), eachdof(rp))
        @test dof == A[I,:]
    end
end

@testset "NormalModeTransformation" begin
    @testset "Transformation matrix, nbeads = $n" for n in [10, 11]
        transform = RingPolymerArrays.NormalModeTransformation{Float64}(n)
        U = transform.U
        @test U'U ≈ I 
        S = Matrix(SymTridiagonal(fill(2, n), fill(-1, n-1)))
        S[end,1] = S[1,end] = -1
        @test U'S*U ≈ Diagonal(4sin.((0:n-1)*π/n) .^2)
    end

    @testset "transform_to/from_normal_modes!" begin
        nbeads = 10
        rp = RingPolymerArray(rand(ndofs, natoms, nbeads))
        initial = deepcopy(rp)
        transform = RingPolymerArrays.NormalModeTransformation{Float64}(nbeads)
        RingPolymerArrays.transform_to_normal_modes!(rp, transform)
        @test !(rp ≈ initial)
        RingPolymerArrays.transform_from_normal_modes!(rp, transform)
        @test rp ≈ initial
    end
end

@testset "get_centroid" begin
    r = rand(3, 4, 6)
    centroid = dropdims(mean(r; dims=3); dims=3)
    @test centroid ≈ RingPolymerArrays.get_centroid(r)

    r = RingPolymerArray(rand(3, 4, 6))
    centroid = dropdims(mean(r; dims=3); dims=3)
    @test centroid ≈ RingPolymerArrays.get_centroid(r)
end

@testset "Operations" begin
    rp = RingPolymerArray(A; nbeads)
    B = zeros(Int, size(A)..., nbeads)
    for i in axes(B,3)
        copy!(view(B,:,:,i), A)
    end
    @test (rp + rp) isa RingPolymerArray
    @test RingPolymerArray(B + B) == (rp + rp)
    @test (B - B) == (rp - rp)
    @test (B .* B) == (rp .* rp)
    @test (B ./ B) == (rp ./ rp)
end

@testset "Operations w/ classical atoms" begin
    A = rand(ndofs, natoms, nbeads)
    rp = RingPolymerArray(A; classical=[1])
    B = copy(A)
    @test (rp + rp) isa RingPolymerArray
    @test (B + B) ≈ (rp + rp)
    @test (B - B) ≈ (rp - rp)
    @test (B .* B) ≈ (rp .* rp)
    @test (B ./ B) == (rp ./ rp)
end
