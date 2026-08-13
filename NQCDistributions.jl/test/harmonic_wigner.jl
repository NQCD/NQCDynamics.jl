using Test
using NQCDistributions

ω = 3
β = 5
m = 10
mom = MomentumHarmonicWigner(ω, β, m)
pos = PositionHarmonicWigner(ω, β, m)
vel = VelocityHarmonicWigner(ω, β, m)

@testset "Correct norm" begin
    norm = tanh(β*ω/2) / π
    @test norm ≈ 1 / (mom.σ * pos.σ * 2π)
    @test norm ≈ 1 / (vel.σ * pos.σ * 2π) / m
end

@testset "Energy comparison" begin
    energy_expectation = ω * (1/(exp(ω*β) - 1) + 1/2)

    velocity_kinetic() = m/2*rand(vel)^2
    momentum_kinetic() = rand(mom)^2/2m

    function wigner_energy_expectation(kinetic, n)
        e = 0.0
        for _=1:n
            e += kinetic() + rand(pos)^2 * m *  ω^2/2
        end
        return e / n
    end

    @test wigner_energy_expectation(velocity_kinetic, 100000) ≈ energy_expectation atol=1e-1
    @test wigner_energy_expectation(momentum_kinetic, 100000) ≈ energy_expectation atol=1e-1
end
