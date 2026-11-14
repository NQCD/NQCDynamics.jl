using Test
using NQCModels
using LinearAlgebra: eigvals

include("test_utils.jl")

innermodel = DoubleWell()

@testset "State index out of bounds" begin
    @test_throws DomainError AdiabaticStateSelector(innermodel, 0)
    @test_throws DomainError AdiabaticStateSelector(innermodel, 3)
end

@testset "Potential, state: $selected_state" for selected_state = 1:2
    model = AdiabaticStateSelector(innermodel, selected_state)
    r = randn(1,1)

    diabatic_potential = potential(innermodel, r)
    correct_value = eigvals(diabatic_potential)[selected_state]
    new_value = potential(model, r)
    @test correct_value ≈ new_value
end

@testset "Potential, state: $selected_state" for selected_state = 1:2
    model = AdiabaticStateSelector(innermodel, selected_state)
    @test test_model(model, 1)
end


