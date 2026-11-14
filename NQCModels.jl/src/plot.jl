
module Plot

using RecipesBase: RecipesBase, @recipe, @series, @userplot
using LinearAlgebra: eigvals, diag, eigvecs

using Unitful: @u_str, uconvert
using UnitfulAtomic: UnitfulAtomic

using ..NQCModels: potential, derivative, nstates, state_independent_potential
using ..ClassicalModels: ClassicalModel
using ..QuantumModels: QuantumModel

@recipe function f(x, model::ClassicalModel)
    V = zeros(size(x))
    D = zeros(size(x))
    for i=1:length(x)
        V[i] = potential(model, hcat(x[i]))
        D[i] = derivative(model, hcat(x[i]))[1]
    end

    xguide --> "r"

    @series begin
        label := "V(r)"
        x .* u"bohr", V .* u"hartree"
    end

    @series begin
        label := "dV(r)dr"
        x .* u"bohr", D .* u"hartree / bohr"
    end
end

@recipe function f(x, model::QuantumModel; adiabats=true, diabats=true, coupling=false, atomic=true, include_diagonal=true)
    eigs = zeros(length(x), nstates(model))
    diabatic = zeros(length(x), nstates(model))
    couplings = zeros(length(x), nstates(model), nstates(model))
    for i=1:length(x)
        V = Matrix(potential(model, hcat(x[i])))
        state_independent = state_independent_potential(model, hcat(x[i]))
        if include_diagonal
            eigs[i,:] .= eigvals(V) .+ state_independent
            diabatic[i,:] .= diag(V) .+ state_independent
        else
            eigs[i,:] .= eigvals(V)
            diabatic[i,:] .= diag(V)
        end
        couplings[i,:,:] .= V
    end

    xguide := "r"
    yguide := "V(r)"

    if adiabats
        for i=1:nstates(model)
            @series begin
                linecolor := :black
                label := i==1 ? "Classical" : ""
                if atomic
                    x .* u"bohr", eigs[:,i] .* u"hartree"
                else
                    uconvert.(u"Å", x .* u"bohr"), uconvert.(u"eV", eigs[:,i] .* u"hartree")
                end
            end
        end
    end

    if diabats
        for i=1:nstates(model)
            @series begin
                linecolor := "#FF1F5B"
                label := i==1 ? "Diabatic" : ""
                if atomic
                    x .* u"bohr", diabatic[:,i] .* u"hartree"
                else
                    uconvert.(u"Å", x .* u"bohr"), uconvert.(u"eV", diabatic[:,i] .* u"hartree")
                end
            end
        end
    end

    if coupling
        for i=1:nstates(model)
            for j=i+1:nstates(model)
                @series begin
                    linecolor := "#009ADE"
                    label := (i == 1 && j == 2) ? "Diabatic coupling" : ""
                    if atomic
                        x .* u"bohr", couplings[:,i,j] .* u"hartree"
                    else
                        uconvert.(u"Å", x .* u"bohr"), uconvert.(u"eV", couplings[:,i,j] .* u"hartree")
                    end
                end
            end
        end
    end
end

@userplot PlotClassicalGradient

@recipe function f(p::PlotClassicalGradient; states=1:1, atomic=true)
    x, model = p.args

    legend := false

    gradient = zeros(length(x), nstates(model), nstates(model))
    for i=1:length(x)
        V = potential(model, hcat(x[i]))
        U = eigvecs(V)
        D = derivative(model, hcat(x[i]))
        gradient[i,:,:] .= U' * D[1] * U
    end

    for i in states
        for j=i:nstates(model)
            @series begin
                linecolor := i == j ? :black : :red
                if atomic
                    x .* u"bohr", (gradient[:,j,i] .* u"hartree") .^2
                else
                    uconvert.(u"Å", x .* u"bohr"), uconvert.(u"eV^2/Å^2", (gradient[:,j,i] .* u"hartree/bohr") .^2)
                end
            end
        end
    end
end

end # module
