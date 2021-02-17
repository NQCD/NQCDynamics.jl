using LinearAlgebra
using RecipesBase

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

    xguide --> "r /a₀"

    @series begin
        label := "V(r) /Eₕ"
        x, potential
    end

    @series begin
        label := "dV(r)dr /Eₕ/a₀"
        x, derivative
    end
end

@recipe function f(x, model::DiabaticModel)
    V = Hermitian(zeros(model.n_states, model.n_states))
    eigs = zeros(length(x), model.n_states)
    diabats = zeros(length(x), model.n_states)
    for i=1:length(x)
        potential!(model, V, hcat(x[i]))
        eigs[i,:] .= eigvals(V)
        diabats[i,:] .= diag(V)
    end

    legend --> false
    xguide --> "r /a₀"
    yguide --> "V(r) /eV"

    @series begin
        linecolor := :black
        x, au_to_eV.(eigs)
    end

    @series begin
        linecolor := :red
        x, au_to_eV.(diabats)
    end
end
