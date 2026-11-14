using Test
using NQCModels
using LinearAlgebra: Hermitian, eigvals, diag

include("test_utils.jl")

# Minimal model that already provides an adiabatic (diagonal) Hamiltonian
struct ToyAdiabaticModel <: NQCModels.QuantumModels.QuantumModel{Adiabatic} end
NQCModels.nstates(::ToyAdiabaticModel) = 2
NQCModels.ndofs(::ToyAdiabaticModel) = 1
function NQCModels.potential!(::ToyAdiabaticModel, V::Hermitian, R::AbstractMatrix)
    fill!(V.data, 0.0)
    V.data[1, 1] = sum(R)
    V.data[2, 2] = 2 * sum(R) + 1
    return nothing
end
function NQCModels.derivative!(
    ::ToyAdiabaticModel,
    D::AbstractMatrix{<:Hermitian},
    R::AbstractMatrix,
)
    for I in eachindex(D)
        fill!(D[I].data, 0.0)
        D[I].data[1, 1] = 1.0
        D[I].data[2, 2] = 2.0
    end
    return nothing
end

@testset "hamiltonian_type" begin
    @test hamiltonian_type(DoubleWell()) == Diabatic()
    @test hamiltonian_type(ErpenbeckThoss(Γ = 0.1)) == Diabatic()
    @test hamiltonian_type(TullyModelOne()) == Diabatic()
    @test hamiltonian_type(
        AndersonHolstein(ErpenbeckThoss(Γ = 0.1), TrapezoidalRule(10, -1, 1)),
    ) == Diabatic()
    @test hamiltonian_type(
        WideBandBath(ErpenbeckThoss(Γ = 0.1); step = 0.1, bandmin = -1.0, bandmax = 1.0),
    ) == Diabatic()
    @test hamiltonian_type(ToyAdiabaticModel()) == Adiabatic()
end

@testset "StateSelector bounds" begin
    @test_throws DomainError StateSelector(DoubleWell(), 0, Adiabatic)
    @test_throws DomainError StateSelector(DoubleWell(), 3, Diabatic)
end

@testset "StateSelector of diabatic model, state $s" for s = 1:2
    inner = DoubleWell()
    r = randn(1, 1)
    V = potential(inner, r)

    adiabatic = StateSelector(inner, s, Adiabatic)
    @test potential(adiabatic, r) ≈ eigvals(V)[s]
    @test test_model(adiabatic, 1)

    diabatic = StateSelector(inner, s, Diabatic)
    @test potential(diabatic, r) ≈ V[s, s]
    @test test_model(diabatic, 1)
end

@testset "StateSelector of adiabatic model" begin
    r = hcat(0.3)
    @test potential(StateSelector(ToyAdiabaticModel(), 1, Adiabatic), r) ≈ 0.3
    @test derivative(StateSelector(ToyAdiabaticModel(), 1, Adiabatic), r) ≈ hcat(1.0)
    # Known issue carried over from NQCD/NQCModels.jl#76: linear indexing V[state] instead of V[state, state]
    @test_broken potential(StateSelector(ToyAdiabaticModel(), 2, Adiabatic), r) ≈ 1.6
    @test_broken derivative(StateSelector(ToyAdiabaticModel(), 2, Adiabatic), r) ≈ hcat(2.0)
end

@testset "ReducedQuantumModel" begin
    inner = ThreeStateMorse()
    r = hcat(3.0)
    V = potential(inner, r)

    @test_throws DomainError ReducedQuantumModel(inner, [0, 1], Diabatic)
    @test_throws ArgumentError ReducedQuantumModel(inner, [1, 1], Diabatic)

    reduced = ReducedQuantumModel(inner, [1, 3], Diabatic)
    @test nstates(reduced) == 2
    @test hamiltonian_type(reduced) == Diabatic()
    @test potential(reduced, r) ≈ V[[1, 3], [1, 3]]
    @test test_model(reduced, 1)

    reduced_adiabatic = ReducedQuantumModel(inner, [1, 2], Adiabatic)
    @test hamiltonian_type(reduced_adiabatic) == Adiabatic()
    @test diag(potential(reduced_adiabatic, r)) ≈ eigvals(V)[1:2]
    @test size(derivative(reduced_adiabatic, r)) == size(r)
end
