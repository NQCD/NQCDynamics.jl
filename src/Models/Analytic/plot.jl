export plot_diabats
export plot_adiabats
export plot_coupling
export plot_eigenvalues

using Plots

"""
    plot_diabats(model::Model, domain::Tuple, points::Integer=75)

Plot the diagonal elements of the potential energy matrix.
"""
function plot_diabats(model::Model, domain::Tuple, points::Integer=75)
    x = range(domain..., length=points)
    V = model.get_potential.(x)
    V0 = model.get_V0.(x)
    n_states = size(V[1])[1]
    p = plot(legend=false)
    for i=1:n_states
        y = [a[i, i] for a in V] .+ V0
        plot!(x, y)
    end
    return p
end

"""
    plot_adiabats(model::Model, domain::Tuple, points::Integer=75)

Plot the eigenvalues of the potential energy matrix.
"""
function plot_adiabats(model::Model, domain::Tuple, points::Integer=75)
    x = range(domain..., length=points)
    V = model.get_potential.(x)
    V0 = model.get_V0.(x)
    vals = eigvals.(V)
    n_states = length(vals[1])
    p = plot(legend=false)
    for i=1:n_states
        y = [a[i] for a in vals] .+ V0
        plot!(x, y)
    end
    return p
end

"""
    plot_coupling(model::Model, domain::Tuple, points::Integer=75)

Plot the off-diagonal elements of the potential energy matrix.
"""
function plot_coupling(model::Model, domain::Tuple, points::Integer=75)
    x = range(domain..., length=points)
    V = model.get_potential.(x)
    n_states = size(V[1])[1]
    p = plot(legend=false)
    for i=1:n_states
        for j=i+1:n_states
            y = [a[i, j] for a in V]
            plot!(x, y, linestyle=:dash)
        end
    end
    return p
end

"""
    plot_eigenvalues(model::Model, R::Real)

Plot the eigenvalues of the potential energy matrix at a single point
"""
function plot_eigenvalues(model::Model, R::Real)
    V = model.get_potential(R)
    V0 = model.get_V0(R)
    plot(V0 .+ eigvals(V), legend=false)
end