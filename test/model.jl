using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Models
using LinearAlgebra

R = rand(3, 10)

@testset "Harmonic" begin
    model = Harmonic()
    model.potential!([0.0], R)
    model.derivative!(zero(R), R)
end

@testset "Free" begin
    model = Free()
    model.potential!([0.0], R)
    model.derivative!(zero(R), R)
end

@testset "DoubleWell" begin
    R = rand(1,1)
    model = DoubleWell()
    model.potential!(Hermitian(zeros(model.n_states, model.n_states)), R)
    model.derivative!([Hermitian(zeros(model.n_states, model.n_states))]', R)
end

@testset "TullyModelOne" begin
    R = rand(1,1)
    model = TullyModelOne()
    model.potential!(Hermitian(zeros(model.n_states, model.n_states)), R)
    model.derivative!([Hermitian(zeros(model.n_states, model.n_states))]', R)
end

@testset "TullyModelTwo" begin
    R = rand(1,1)
    model = TullyModelTwo()
    model.potential!(Hermitian(zeros(model.n_states, model.n_states)), R)
    model.derivative!([Hermitian(zeros(model.n_states, model.n_states))]', R)
end

@testset "PdH" begin
    model = PdH([:Pd, :Pd, :H], PeriodicCell([10 0 0; 0 10 0; 0 0 10]), 10.0)
    R = rand(3, 3) * 10
    V = [0.0]
    D = zero(R)
    model.potential!(V, R)
    model.derivative!(D, R)
    h = 1e-4
    V1 = [0.0]
    for i=1:length(R)
        R[i] += h
        model.potential!(V1, R)
        @test D[i] ≈ (V1 - V)[1] / h rtol=1e-1
        R[i] -= h
    end
end