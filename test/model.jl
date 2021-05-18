using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Models
using LinearAlgebra
using FiniteDiff

function finite_difference_gradient(model::AdiabaticModel, R)
    f(x) = potential(model, x)[1]
    FiniteDiff.finite_difference_gradient(f, R)
end

function finite_difference_gradient(model::DiabaticModel, R)
    f(x, i, j) = potential(model, x)[i,j]
    grad = [Hermitian(zeros(model.n_states, model.n_states)) for i in CartesianIndices(R)]
    for k in eachindex(R)
        for i=1:model.n_states
            for j=1:model.n_states
                grad[k].data[i,j] = FiniteDiff.finite_difference_gradient(x->f(x,i,j), R)[k]
            end
        end
    end
    grad
end

function test_model(model::AdiabaticModel, DoFs, atoms)
    R = rand(DoFs, atoms)
    D = derivative(model, R)
    @test finite_difference_gradient(model, R) ≈ D
end

function test_model(model::DiabaticModel, DoFs, atoms)
    R = rand(DoFs, atoms)
    D = derivative(model, R)
    @test finite_difference_gradient(model, R) ≈ D rtol=1e-3
end

function test_model(model::AdiabaticFrictionModel, DoFs, atoms)
    R = rand(DoFs, atoms)
    D = derivative(model, R)
    friction(model, R)
    @test finite_difference_gradient(model, R) ≈ D
end

@testset "Harmonic" begin
    model = Harmonic()
    test_model(model, 3, 10)
end

@testset "DiatomicHarmonic" begin
    model = DiatomicHarmonic()
    test_model(model, 3, 2)

    R = [0 0; 0 0; 1 0]
    @test Models.energy(model, R) ≈ 0
    R = [sqrt(3) 0; sqrt(3) 0; sqrt(3) 0]
    @test Models.energy(model, R) ≈ 2
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

@testset "ThreeStateMorse" begin
    model = ThreeStateMorse()
    test_model(model, 1, 1)
end

@testset "DebyeSpinBoson" begin
    model = DebyeSpinBoson(10)
    test_model(model, 1, 10)
end

@testset "DebyeBosonBath" begin
    model = DebyeBosonBath(10)
    test_model(model, 1, 10)
end

@testset "JuLIP" begin
    import JuLIP
    atoms = NonadiabaticMolecularDynamics.Atoms{Float64}([:H, :H])
    x = [10.0 0 0]
    y = [0 10.0 0]
    z = [0 0 10.0]
    model = JuLIPModel(atoms, PeriodicCell(vcat(x,y,z)), JuLIP.StillingerWeber())
    test_model(model, 3, 2)
end

# Commented because it requires a cube file to test.
# @testset "LDFA" begin
#     import PyCall
#     free = Free()
#     atoms = Atoms([:H, :H])
#     cube_file_here = "This should be a cube file"
#     model = Models.LDFA.CubeLDFAModel(free, cube_file_here, atoms, [2])
#     test_model(model, 3, 2)
# end
