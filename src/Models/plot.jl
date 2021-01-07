using Plots
using LinearAlgebra

function Plots.plot(x, model::AdiabaticModel)
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
    plot(x, potential, label="potential!")
    plot!(x, derivative, label="derivative!")
end

function Plots.plot(x, model::DiabaticModel)
    V = Hermitian(zeros(model.n_states, model.n_states))
    eigs = zeros(length(x), model.n_states)
    diabats = zeros(length(x), model.n_states)
    for i=1:length(x)
        potential!(model, V, hcat(x[i]))
        eigs[i,:] .= eigvals(V)
        diabats[i,:] .= diag(V)
    end
    plot(x, eigs, color="black", legend=false)
    plot!(x, diabats, color="red")
end