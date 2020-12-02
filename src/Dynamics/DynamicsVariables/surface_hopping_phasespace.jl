using RecursiveArrayTools

export SurfaceHoppingPhasespace
export get_density_matrix
export get_adiabatic_state

"""
    SurfaceHoppingPhasespace{T} <: DynamicalVariables{T}
    
Phasespace that includes a density matrix and a discrete state variable.
"""
struct SurfaceHoppingPhasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}}
    state::UInt
end

"""
    SurfaceHoppingPhasespace(R::Vector{T}, P::Vector{T}, n_states::Integer, state::Integer) where {T<:AbstractFloat}
    
Constructor for nonequilibrium surface hopping phasespace initialised on adiabatic state `state`.
"""
function SurfaceHoppingPhasespace(R::Vector{T}, P::Vector{T}, n_states::Integer, state::Integer) where {T<:AbstractFloat}
    z = Phasespace(R, P) 
    σ = zeros(Complex{T}, n_states, n_states)
    σ[state, state] = 1

    SurfaceHoppingPhasespace{T}(ArrayPartition(z.x, σ), state)
end

get_positions(z::SurfaceHoppingPhasespace) = @view z.x.x[1][:,1]
get_momenta(z::SurfaceHoppingPhasespace) = @view z.x.x[1][:,2]
get_density_matrix(z::SurfaceHoppingPhasespace) = z.x.x[2]
get_adiabatic_state(z::SurfaceHoppingPhasespace) = z.state