
"""
    ElectronicDistribution{S}

Abstract type for distributions of electronic degrees of freedom only.
"""
abstract type ElectronicDistribution{S} end


"""
    PureState{S} <: ElectronicDistribution{S}

Electronic distribution for representing a system confined to a single state.
"""
struct PureState{S} <: ElectronicDistribution{S}
    state::Int
    statetype::S
end
PureState(state) = PureState(state, Diabatic())

function density_matrix(d::PureState, nstates)
    density = zeros(nstates, nstates)
    density[d.state, d.state] = 1
    return density
end

"""
    MixedState{T,S} <: ElectronicDistribution{S}

Electronic distribution for representing a mixed state with non-zero population in multiple states.
"""
struct MixedState{T,S} <: ElectronicDistribution{S}
    populations::Vector{T}
    statetype::S
end
MixedState(state) = MixedState(state, Diabatic())

function density_matrix(d::MixedState)
    density = zeros(length(d.populations), length(d.populations))
    density[diagind(density)] .= d.populations
    return density
end

"""
    FermiDiracState{S,T,A} <: ElectronicDistribution{S}

Electronic distribution for Fermions following Fermi-Dirac distribution.
"""
struct FermiDiracState{S,T,A} <: ElectronicDistribution{S}
    fermi_level::T
    β::T
    statetype::S
    available_states::A
end
function FermiDiracState(fermi_level, temperature; statetype=Adiabatic(), available_states=Colon())
    return FermiDiracState(austrip(fermi_level), 1/austrip(temperature), statetype, available_states)
end
fermi(ϵ, μ, β) = 1 / (1 + exp(β*(ϵ - μ)))

function density_matrix(d::FermiDiracState, eigenvalues)
    density = zeros(length(eigenvalues), length(eigenvalues))
    density[diagind(density)] .= fermi.(eigenvalues, d.fermi_level, d.β)
    return density
end

# ------------------------------- New Electronic State Struct/Type ------------------------------- #
"""
    NonEqState{S,T,A} <: ElectronicDistribution{S}

Electronic distribution for Fermions following a supplied non-equilibrium distribution.
"""
struct NonEqState{S,T,A} <: ElectronicDistribution{S}
    dis_spline::T # non-equilibrium distribution spline function
    dos_spline::T # density of states distribution spline function
    statetype::S # usually the `adiabatic()` or `diabatic()` labels
    available_states::A # need this, in general will default to `Colon()` as with above.
end
function NonEqState(dis_spline::LinearInterpolation, dos_spline::LinearInterpolation; statetype=Adiabatic(), available_states=Colon()) 
    return NonEqState(dis_spline, dos_spline, statetype, available_states)
end
# ------------------------------------------------------------------------------------------------ #



# ---------------------------------------- Density Matrix ---------------------------------------- #



"Singleton type for labelling the density matrix occupations as statistical."
struct Statistical end
"Singleton type for labelling the density matrix occupations as binary."
struct Binary end


"""
    DensityMatrix{T,S} <: ElectronicDistribution{S}

Electronic distribution that requests a density matrix representation of the electronic subsystem.

`diagonal` holds the distribution used to populate the diagonal of the density matrix
(a [`FermiDiracState`](@ref) or a [`NonEqState`](@ref)).
`systemtype` labels the occupations as [`Binary`](@ref) or [`Statistical`](@ref).

    DensityMatrix(fermi_level::Float64, temperature::Float64; statetype=Adiabatic(), available_states=Colon(), systemtype=Binary())
    DensityMatrix(dis_spline::LinearInterpolation, dos_spline::LinearInterpolation; statetype=Adiabatic(), available_states=Colon(), systemtype=Binary())

The first form builds a Fermi-Dirac diagonal (atomic units), the second a non-equilibrium diagonal.
"""
struct DensityMatrix{T,S} <: ElectronicDistribution{S}
    diagonal::T
    systemtype::S
end
function DensityMatrix(fermi_level::Float64, temperature::Float64; statetype=Adiabatic(), available_states=Colon(), systemtype=Binary())
    dist = FermiDiracState(austrip(fermi_level), 1/austrip(temperature), statetype, available_states)
    return DensityMatrix(dist, systemtype)
end
function DensityMatrix(dis_spline::LinearInterpolation, dos_spline::LinearInterpolation; statetype=Adiabatic(), available_states=Colon(), systemtype=Binary())
    dist = NonEqState(dis_spline, dos_spline, statetype, available_states)
    return DensityMatrix(dist, systemtype)
end
