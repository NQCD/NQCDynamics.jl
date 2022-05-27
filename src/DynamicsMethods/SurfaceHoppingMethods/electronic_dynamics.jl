
using SciMLBase: ODEProblem, DiffEqArrayOperator
using LinearAlgebra: rmul!

struct ElectronicParameters{T}
    "Product of velocity and nonadiabatic coupling"
    dynamical_coupling::Matrix{Complex{T}} 
    "Eigenvalues used to propagate"
    eigenvalues::Vector{Complex{T}}
end

function set_eigenvalues!(p::ElectronicParameters, eigenvalues::AbstractVector)
    copy!(p.eigenvalues, eigenvalues)
end

function set_dynamical_coupling!(p::ElectronicParameters, nonadiabatic_coupling::AbstractArray, velocity::AbstractArray)
    fill!(p.dynamical_coupling, zero(eltype(p.dynamical_coupling)))
    for I in eachindex(nonadiabatic_coupling, velocity)
        p.dynamical_coupling .+= nonadiabatic_coupling[I] .* velocity[I]
    end
end

function set_dynamical_coupling!(p::ElectronicParameters, dynamical_coupling::AbstractArray)
    copy!(p.dynamical_coupling, dynamical_coupling)
end

function update_operator!(L::DiffEqArrayOperator, p::ElectronicParameters)
    A = L.A
    copy!(A, p.dynamical_coupling)
    rmul!(A, -1)
    for (i, I) in zip(eachindex(p.eigenvalues), diagind(A))
        A[I] = -im * p.eigenvalues[i]
    end
end

function ElectronicODEProblem(c0, tspan, nstates::Integer)
    p = ElectronicParameters(zeros(eltype(c0), nstates, nstates), zeros(eltype(c0), nstates))
    A = DiffEqArrayOperator(zeros(ComplexF64, nstates, nstates))
    return ODEProblem(A, c0, tspan, p)
end
