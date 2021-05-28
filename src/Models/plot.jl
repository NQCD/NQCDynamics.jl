using LinearAlgebra
using RecipesBase
using Unitful, UnitfulRecipes

@recipe function f(x, model::AdiabaticModel)
    V = [0.0]
    D = hcat(0.0)
    potential = zeros(size(x))
    derivative = zeros(size(x))
    for i=1:length(x)
        potential!(model, V, hcat(x[i]))
        derivative!(model, D, hcat(x[i]))
        potential[i] = V[1]
        derivative[i] = D[1]
    end

    xguide --> "r"

    @series begin
        label := "V(r)"
        x .* u"bohr", potential .* u"hartree"
    end

    @series begin
        label := "dV(r)dr"
        x .* u"bohr", derivative .* u"hartree / bohr"
    end
end

@recipe function f(x, model::DiabaticModel; adiabats=true, diabats=true)
    V = Hermitian(zeros(model.n_states, model.n_states))
    eigs = zeros(length(x), model.n_states)
    diabatic = zeros(length(x), model.n_states)
    for i=1:length(x)
        potential!(model, V, hcat(x[i]))
        eigs[i,:] .= eigvals(V)
        diabatic[i,:] .= diag(V)
    end

    legend --> false
    xguide --> "r"
    yguide --> "V(r)"

    if adiabats
        @series begin
            linecolor := :black
            x .* u"bohr", eigs .* u"hartree"
        end
    end

    if diabats
        @series begin
            linecolor := :red
            x .* u"bohr", diabatic .* u"hartree"
        end
    end
end
