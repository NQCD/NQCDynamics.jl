using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Models
using LinearAlgebra
using FiniteDiff

function finite_difference_gradient(model::Union{AdiabaticModel, FrictionModel}, R)
    function f(x)
        out = [0.0]
        potential!(model, out, x)
        out[1]
    end
    FiniteDiff.finite_difference_gradient(f, R)
end

function finite_difference_gradient(model::DiabaticModel, R)
    function f(x, i, j)
        out = Hermitian(zeros(model.n_states, model.n_states))
        potential!(model, out, x)
        out[i,j]
    end
    grad = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:size(R)[1], j=1:size(R)[2]]
    for i=1:model.n_states
        for j=1:model.n_states
            grad[1].data[i,j] = FiniteDiff.finite_difference_gradient(x->f(x,i,j), R)[1]
        end
    end
    grad
end

function test_model(model::AdiabaticModel, DoFs, atoms)
    R = rand(DoFs, atoms)
    V = zeros(1)
    D = zero(R)
    potential!(model, V, R)
    derivative!(model, D, R)
    @test finite_difference_gradient(model, R) ≈ D
end

function test_model(model::DiabaticModel, DoFs, atoms)
    R = rand(DoFs, atoms)
    V = Hermitian(zeros(model.n_states, model.n_states))
    D = [Hermitian(zeros(model.n_states, model.n_states))]'
    potential!(model, V, R)
    derivative!(model, D, R)
    @test finite_difference_gradient(model, R) ≈ D rtol=1e-3
end

function test_model(model::FrictionModel, DoFs, atoms)
    R = rand(DoFs, atoms)
    V = zeros(1)
    D = zero(R)
    F = zeros(DoFs*atoms, DoFs*atoms)
    potential!(model, V, R)
    derivative!(model, D, R)
    @test finite_difference_gradient(model, R) ≈ D
    friction!(model, F, R)
end

@testset "Harmonic" begin
    model = Harmonic()
    test_model(model, 3, 10)
end

@testset "DiatomicHarmonic" begin
    model = DiatomicHarmonic()
    test_model(model, 3, 2)
end

@testset "Free" begin
    model = Free()
    test_model(model, 3, 10)
end

@testset "DoubleWell" begin
    model = DoubleWell()
    test_model(model, 1, 1)
end

@testset "TullyModelOne" begin
    model = TullyModelOne()
    test_model(model, 1, 1)
end

@testset "TullyModelTwo" begin
    model = TullyModelTwo()
    test_model(model, 1, 1)
end

@testset "FrictionHarmonic" begin
    model = FrictionHarmonic()
    test_model(model, 1, 3)
end

@testset "ScatteringAndersonHolstein" begin
    model = ScatteringAndersonHolstein()
    test_model(model, 1, 1)
end

@testset "Scattering1D" begin
    model = Scattering1D()
    test_model(model, 1, 1)
end

@testset "PdH" begin
    model = PdH([:Pd, :Pd, :H], PeriodicCell([1 0 0; 0 1 0; 0 0 1]), 1.0)
    test_model(model, 3, 3)
end
