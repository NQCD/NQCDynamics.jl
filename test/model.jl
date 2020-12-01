using Test
using NonadiabaticMolecularDynamics.Models
using NonadiabaticMolecularDynamics.Models.Analytic

function call_functions(model::Model, R::AbstractFloat)
    model.get_V0(R)
    model.get_D0(R)
    model.get_potential(R)
    model.get_derivative(R)
end


@testset "Create models" begin
    R = 10.0

    model = Harmonic(1, 1, 0)
    call_functions(model, R)

    model = MorsePotential()
    call_functions(model, R)

    model = TullyModelOne()
    call_functions(model, R)

    model = TullyModelTwo()
    call_functions(model, R)

    model = Free()
    call_functions(model, R)

    model = MetallicChain()
    call_functions(model, R)

    model = Eckart(1, 1)
    call_functions(model, R)

    model = DoubleWell()
    call_functions(model, R)

    model = DoubleWell(1, 1, 1, 1)
    call_functions(model, R)

    model = DiffusionMetal()
    call_functions(model, R)
    
    model = ScatteringMetal()
    call_functions(model, R)

end