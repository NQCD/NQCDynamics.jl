using Test
using NQCDynamics
using LinearAlgebra

correction = SurfaceHoppingMethods.DecoherenceCorrectionEDC()

ψ = normalize(rand(ComplexF64, 4))
ψ_original = copy(ψ)

occupied_state = 2
dt = 0.1
E = [0.1, 0.2, 0.3, 0.4]
Ekin = 0.25

@testset "apply_decoherence_correction!" begin
    SurfaceHoppingMethods.apply_decoherence_correction!(
    ψ, correction, occupied_state, dt, E, Ekin
    )

    for i in eachindex(ψ)
        if i != occupied_state
            @test abs2(ψ[i]) < abs2(ψ_original[i])
        else
            @test abs2(ψ[i]) > abs2(ψ_original[i])
        end
    end
    @test norm(ψ_original) ≈ 1
    @test norm(ψ) ≈ 1
end
