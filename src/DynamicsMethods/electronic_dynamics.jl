
using SciMLBase: ODEProblem, DiffEqArrayOperator
using LinearAlgebra: rmul!, diagind
using NQCDynamics: DynamicsUtils

struct ElectronicParameters{T}
    "Product of velocity and nonadiabatic coupling"
    dynamical_coupling::Matrix{Complex{T}} 
    "Eigenvalues used to propagate"
    eigenvalues::Vector{Complex{T}}
    t::Base.RefValue{T}
end

mutable struct DoubleBuffer{A}
    current::A
    next::A
end

function swap!(buffer::DoubleBuffer)
    tmp = buffer.current
    setfield!(buffer, :current, buffer.next)
    setfield!(buffer, :next, tmp)
end

function Base.setproperty!(buffer::DoubleBuffer, name::Symbol, x)
    setfield!(buffer.next, name, x)
end

function Base.getproperty(buffer::DoubleBuffer, name::Symbol)
    if name === :current || name === :next
        return Base.getfield(buffer, name)
    else
        return Base.getfield(buffer.next, name)
    end
end

function update_parameters!(buffer::DoubleBuffer,
    eigenvalues::AbstractVector,
    nonadiabatic_coupling::AbstractMatrix{<:AbstractMatrix},
    velocity::AbstractMatrix,
    t::Real
)
    swap!(buffer)
    copy!(buffer.eigenvalues, eigenvalues)
    fill!(buffer.dynamical_coupling, zero(eltype(buffer.dynamical_coupling)))
    for I in eachindex(nonadiabatic_coupling, velocity)
        for J in eachindex(buffer.dynamical_coupling, nonadiabatic_coupling[I])
            buffer.dynamical_coupling[J] += nonadiabatic_coupling[I][J] * velocity[I]
        end
    end
    buffer.t[] = t
end

function interpolate_eigenvalues!(output::AbstractVector, buffer::DoubleBuffer, t)
    E₀ = buffer.current.eigenvalues
    E₁ = buffer.next.eigenvalues
    t₀ = buffer.current.t[]
    t₁ = buffer.next.t[]
    
    interval_location = (t - t₀) / (t₁ - t₀)
    isnan(interval_location) && (interval_location = 0)
    for i in eachindex(E₀, E₁)
        output[i] = E₀[i] + (E₁[i]-E₀[i]) * interval_location
    end
end

function interpolate_dynamical_coupling!(output::AbstractMatrix, buffer::DoubleBuffer, t)
    d₀ = buffer.current.dynamical_coupling
    d₁ = buffer.next.dynamical_coupling
    t₀ = buffer.current.t[]
    t₁ = buffer.next.t[]
    
    interval_location = (t - t₀) / (t₁ - t₀)
    isnan(interval_location) && (interval_location = 0)
    for i in eachindex(d₀, d₁)
        output[i] = d₀[i] + (d₁[i]-d₀[i]) * interval_location
    end
end

function ElectronicODEProblem(c0, tspan, nstates::Integer)
    p = ElectronicParameters(zeros(eltype(c0), nstates, nstates), zeros(eltype(c0), nstates), Ref(tspan[1]))
    buffer = DoubleBuffer(p, deepcopy(p))

    function update_func(A, u, p, t)
        interpolate_dynamical_coupling!(A, p, t)
        diagonal = @view A[diagind(A)]
        interpolate_eigenvalues!(diagonal, p, t)
        rmul!(A, -1)
        rmul!(diagonal, im)
    end

    A = DiffEqArrayOperator(zeros(ComplexF64, nstates, nstates); update_func)
    return ODEProblem(A, c0, tspan, buffer)
end

struct DensityMatrixPropagationParameters{T}
    buffer::DoubleBuffer{ElectronicParameters{T}}
    propagator::Matrix{Complex{T}}
end

function DensityMatrixODEProblem(u0::AbstractMatrix, tspan, nstates::Integer)

    function f!(du, u, p, t)
        A = p.propagator
        interpolate_dynamical_coupling!(A, p.buffer, t)
        rmul!(A, -1im)

        diagonal = @view A[diagind(A)]
        interpolate_eigenvalues!(diagonal, p.buffer, t)

        DynamicsUtils.commutator!(du, A, u)
        rmul!(du, -im)
    end

    electronic_parameters = ElectronicParameters(
        zeros(eltype(u0), nstates, nstates),
        zeros(eltype(u0), nstates),
        Ref(tspan[1])
    )

    p = DensityMatrixPropagationParameters(
        DoubleBuffer(electronic_parameters, deepcopy(electronic_parameters)),
        zeros(eltype(u0), nstates, nstates)
    )

    return ODEProblem(f!, u0, tspan, p)
end

function update_parameters!(
    p::DensityMatrixPropagationParameters,
    eigenvalues::AbstractVector,
    nonadiabatic_coupling::AbstractMatrix{<:AbstractMatrix},
    velocity::AbstractMatrix,
    t::Real
)
    update_parameters!(p.buffer, eigenvalues, nonadiabatic_coupling, velocity, t)
end
