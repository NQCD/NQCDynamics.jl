using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Models

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