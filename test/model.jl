using Test
using NonadiabaticMolecularDynamics
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
    
    model = PdH([:Pd, :Pd, :H], Atoms.PeriodicCell([10 0 0; 0 10 0; 0 0 10]), 10.0)
    R = rand(Float64, 3, 3) * 10
    V = model.get_V0(R)
    D = model.get_D0(R)
    h = 1e-4
    for i=1:length(R)
        R[i] += h
        V1 = model.get_V0(R)
        @test D[i] ≈ (V1 - V) / h rtol=1e-1
        R[i] -= h
    end

end