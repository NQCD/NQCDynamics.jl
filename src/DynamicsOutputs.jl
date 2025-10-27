
"""
    DynamicsOutputs

Defines a set of functions that can be used to calculate outputs for dynamics simulations.
"""
module DynamicsOutputs

using RingPolymerArrays: get_centroid
using Unitful
using UnitfulAtomic
using LinearAlgebra: norm
using ComponentArrays: ComponentVector

using NQCDynamics: Estimators, DynamicsUtils, Analysis, ndofs

using NQCModels: NQCModels

using .DynamicsUtils: get_positions, get_velocities, get_quantum_subsystem

using ..InitialConditions: QuantisedDiatomic

using NQCBase
using Statistics

"""
    OutputVelocity(sol, i)

Output the velocity at each timestep during the trajectory.
"""
OutputVelocity(sol, i) = [copy(get_velocities(u)) for u in sol.u]
export OutputVelocity

"""
    OutputCentroidVelocity(sol, i)

Output the velocity of the ring polymer centroid at each timestep during the trajectory.
"""
OutputCentroidVelocity(sol, i) = [get_centroid(get_velocities(u)) for u in sol.u]
export OutputCentroidVelocity

function OutputCentroidKineticEnergy(sol, i)
    output = zeros(length(sol.u))
    for i in eachindex(output)
        u = sol.u[i]
        centroid = get_centroid(get_velocities(u))
        output[i] = DynamicsUtils.classical_kinetic_energy(sol.prob.p, centroid)
    end
    return output
end
export OutputCentroidKineticEnergy

"""
    OutputPosition(sol, i)

Output the position at each timestep during the trajectory.
"""
OutputPosition(sol, i) = [copy(get_positions(u)) for u in sol.u]
export OutputPosition

"""
    OutputCentroidPosition(sol, i)

Output the position of the ring polymer centroid at each timestep during the trajectory.
"""
OutputCentroidPosition(sol, i) = [get_centroid(get_positions(u)) for u in sol.u]
export OutputCentroidPosition

"""
    OutputPotentialEnergy(sol, i)

Output the adiabatic potential energy at each timestep during the trajectory.
"""
OutputPotentialEnergy(sol, i) = DynamicsUtils.classical_potential_energy.(sol.prob.p, sol.u)
export OutputPotentialEnergy

"""
    OutputTotalEnergy(sol, i)

Evaluate the classical Hamiltonian at each timestep during the trajectory.
"""
OutputTotalEnergy(sol, i) = DynamicsUtils.classical_hamiltonian.(sol.prob.p, sol.u)
export OutputTotalEnergy

"""
    OutputKineticEnergy(sol, i)

Evaluate the classical kinetic energy at each timestep during the trajectory.
"""
OutputKineticEnergy(sol, i) = DynamicsUtils.classical_kinetic_energy.(sol.prob.p, sol.u)
export OutputKineticEnergy

OutputFinalKineticEnergy(sol, i) =
    DynamicsUtils.classical_kinetic_energy(sol.prob.p, last(sol.u))
export OutputFinalKineticEnergy

"""
    OutputSubsetKineticEnergy(sol, i)

Evaluate the classical kinetic energy of a subset of the entire system at each save step.

The subset is defined by `OutputSubsetKineticEnergy(indices)`.
"""
struct OutputSubsetKineticEnergy{T}
    indices::T
end
function (output::OutputSubsetKineticEnergy)(sol, i)
    return map(
        x -> DynamicsUtils.classical_kinetic_energy(
            sol.prob.p.atoms.masses[output.indices],
            x[:, output.indices],
        ),
        [DynamicsUtils.get_velocities(i) for i in sol.u],
    )
end
export OutputSubsetKineticEnergy

"""
Evaluate the classical kinetic energy of a subset of the entire system at the end of the simulation.

The subset is defined by `OutputSubsetKineticEnergy(indices)`.
"""
struct OutputFinalSubsetKineticEnergy
    indices::Any
end
function (output::OutputFinalSubsetKineticEnergy)(sol, i)
    return DynamicsUtils.classical_kinetic_energy(
        sol.prob.p.atoms.masses[output.indices],
        DynamicsUtils.get_velocities(last(sol.u))[:, output.indices],
    )
end
export OutputFinalSubsetKineticEnergy

OutputSpringEnergy(sol, i) = DynamicsUtils.classical_spring_energy.(sol.prob.p, sol.u)
export OutputSpringEnergy

"""
    OutputDynamicsVariables(sol, i)

Output all of the dynamics variables at each timestep during the trajectory.
"""
OutputDynamicsVariables(sol, i) = copy.(sol.u)
export OutputDynamicsVariables

"""
    OutputQuantumSubsystem(sol, i)

Output the quantum subsystem at each timestep during the trajectory.
Usually this will refer to a wavefunction or density matrix but will depend on the particular dynamics method.
"""
OutputQuantumSubsystem(sol, i) = [copy(get_quantum_subsystem(u)) for u in sol.u]
export OutputQuantumSubsystem

"""
    OutputMappingPosition(sol, i)

Output the position mapping variables at each timestep during the trajectory.
"""
OutputMappingPosition(sol, i) =
    [copy(DynamicsUtils.get_mapping_positions(u)) for u in sol.u]
export OutputMappingPosition

"""
    OutputMappingMomentum(sol, i)

Output the momentum mapping variable at each timestep during the trajectory.
"""
OutputMappingMomentum(sol, i) = [copy(DynamicsUtils.get_mapping_momenta(u)) for u in sol.u]
export OutputMappingMomentum

"""
    OutputDiscreteState(sol, i)

Output the discrete state variable at each timestep during the trajectory.
This is used for surface hopping simulations and returns the variable that determines the currently occupied adiabatic state.

Requires that the dynamics variable has a field `state`.

Use [`OutputDiabaticPopulation`](@ref) or [`OutputAdiabaticPopulation`](@ref) to get the population estimators.
"""
OutputDiscreteState(sol, i) = [copy(u.state) for u in sol.u]
export OutputDiscreteState

"""
    OutputDiabaticPopulation(sol, i)

Output the diabatic population at each timestep during the trajectory.
"""
OutputDiabaticPopulation(sol, i) = Estimators.diabatic_population.(sol.prob.p, sol.u)
export OutputDiabaticPopulation

"""
    OutputTotalDiabaticPopulation(sol, i)

Output the total diabatic population at eah timestep during the trajectory.
"""
OutputTotalDiabaticPopulation(sol, i) =
    sum.(Estimators.diabatic_population(sol.prob.p, sol.u))
export OutputTotalDiabaticPopulation

"""
    OutputAdiabaticPopulation(sol, i)

Output the adiabatic population at each timestep during the trajectory.
"""
OutputAdiabaticPopulation(sol, i) = Estimators.adiabatic_population.(sol.prob.p, sol.u)
export OutputAdiabaticPopulation

"""
    OutputTotalAdiabaticPopulation(sol, i)

Output the total adiabatic population at each timestep during the trajectory.
"""
OutputTotalAdiabaticPopulation(sol, i) =
    sum.(Estimators.adiabatic_population.(sol.prob.p, sol.u))
export OutputTotalAdiabaticPopulation

"""
Output the first point of each trajectory in DynamicsVariables format. (Useful when using distributions for initial conditions.)
"""
OutputInitial(sol, i) = first(sol.u)
export OutputInitial

"""
Output the end point of each trajectory.
"""
OutputFinal(sol, i) = last(sol.u)
export OutputFinal

"""
Outputs the final time point of a trajectory. This is useful if simulations are terminated by a callback.
"""
OutputFinalTime(sol, i) = last(sol.t)
export OutputFinalTime


"""
Output the total number of surface hops during the trajectory
"""
function OutputSurfaceHops(sol, i)::Int
    nhops = 0
    for i = 1:length(sol.u)-1
        if sol.u[i].state != sol.u[i+1].state
            nhops += 1
        end
    end
    return nhops
end
export OutputSurfaceHops

"""
Output a 1 if the molecule has dissociated, 0 otherwise.
"""
struct OutputDissociation{T}
    "The maximum distance at which the two atoms can be considered bonded."
    distance::T
    "The indices of the two atoms in the molecule of interest."
    atom_indices::Tuple{Int,Int}
    OutputDissociation(distance, atom_indices) =
        new{typeof(distance)}(austrip(distance), atom_indices)
end
export OutputDissociation

function (output::OutputDissociation)(sol, i)
    R = DynamicsUtils.get_positions(last(sol.u))
    dissociated =
        norm(R[:, output.atom_indices[1]] .- R[:, output.atom_indices[2]]) > output.distance
    return dissociated ? 1 : 0
end

"""
Outputs a 1 if the molecule has moved a certain distance from all atoms in the surface, 0 otherwise.

Distance is defined by the projected distance from the highest non-adsorbate atom along the surface normal.
"""
struct OutputSurfaceDesorption
    distance::Number
    adsorbate_indices::Vector{Int}
    surface_normal::AbstractVector
    OutputSurfaceDesorption(distance, adsorbate_indices; surface_normal = [0, 0, 1]) =
        new(distance, adsorbate_indices, surface_normal)
end

function (osd::OutputSurfaceDesorption)(sol, i)
    desorption_frame = findfirst([
        Analysis.Diatomic.surface_distance_condition(
            i,
            osd.adsorbate_indices,
            sol.prob.p;
            surface_distance_threshold = osd.distance,
            surface_normal = osd.surface_normal,
        ) for i in sol.u
    ])
    return desorption_frame === nothing ? 0 : 1
end

export OutputSurfaceDesorption

"""
Output the vibrational and rotational quantum numbers of the final image.
"""
struct OutputQuantisedDiatomic{H,V}
    height::H
    normal_vector::V
end
OutputQuantisedDiatomic(; height = 10, normal_vector = [0, 0, 1]) =
    OutputQuantisedDiatomic(height, normal_vector)
export OutputQuantisedDiatomic

function (output::OutputQuantisedDiatomic)(sol, i)
    final = last(sol.u)
    ν, J = QuantisedDiatomic.quantise_diatomic(
        sol.prob.p,
        DynamicsUtils.get_velocities(final),
        DynamicsUtils.get_positions(final);
        height = output.height,
        normal_vector = output.normal_vector,
    )
    return (ν, J)
end

"""
Output a `ComponentVector` with fields `reflection` and `transmission` containing
the probability of the outcome.
Each index in the arrays refers to the adiabatic state.
"""
struct OutputStateResolvedScattering1D{S}
    sim::S
    type::Symbol
end
function (output::OutputStateResolvedScattering1D)(sol, i)
    final = last(sol.u) # get final configuration from trajectory
    if output.type == :adiabatic
        populations = Estimators.adiabatic_population(output.sim, final)
    elseif output.type == :diabatic
        populations = Estimators.diabatic_population(output.sim, final)
    else
        throw(ArgumentError("$(output.type) not recognised.
            Only `:diabatic` or `:adiabatic` accepted."))
    end
    output = ComponentVector(
        reflection = zeros(NQCModels.nstates(output.sim)),
        transmission = zeros(NQCModels.nstates(output.sim)),
    )
    x = DynamicsUtils.get_positions(final)[1]
    if x > 0 # If final position past 0 then we count as transmission
        output.transmission .= populations
    else # If final position left of 0 then we count as reflection
        output.reflection .= populations
    end
    return output
end
export OutputStateResolvedScattering1D

"""
Outputs the desorption angle in degrees (relative to the surface normal) if a desorption event is detected.
"""
struct OutputDesorptionAngle{I<:Vector{Int},S<:Vector{Float64},D}
    indices::I
    surface_normal::S
    surface_distance_threshold::D
end
"""
    OutputDesorptionAngle(indices; surface_normal=[0, 0, 1], surface_distance_threshold=austrip(5.0u"Å"))

Outputs the desorption angle in degrees (relative to the surface normal) if a desorption event is detected.
Use `surface_normal` to define the direction "away" from the surface. Most commonly, this would be in positive z direction.

A desorption is detected if the centre of mass of the molecule defined with `indices` is above `surface_distance_threshold` from the closest surface atom.
This is calculated with respect to `surface_normal` and will take into account periodic boundary conditions.
"""
OutputDesorptionAngle(
    indices;
    surface_normal = [0, 0, 1],
    surface_distance_threshold = austrip(5.0u"Å"),
) = OutputDesorptionAngle(
    indices,
    convert(Vector{Float64}, surface_normal),
    surface_distance_threshold,
)
export OutputDesorptionAngle

"""
    (output::OutputDesorptionAngle)(sol, i)

Outputs the desorption angle in degrees (relative to the surface normal) if a desorption event was detected.
"""
function (output::OutputDesorptionAngle)(sol, i)
    return Analysis.Diatomic.get_desorption_angle(
        sol.u,
        output.indices,
        sol.prob.p;
        surface_normal = output.surface_normal,
        surface_distance_threshold = output.surface_distance_threshold,
    )
end

struct OutputDesorptionTrajectory{I<:Vector{Int},N<:Vector{Float64},D,F<:Int}
    indices::I
    surface_normal::N
    surface_distance_threshold::D
    extra_frames::F
end
"""
    OutputDesorptionTrajectory(indices; surface_normal=[0, 0, 1], surface_distance_threshold=austrip(5.0u"Å"), extra_frames=0)

Like OutputDynamicsVariables, but only saves parts of the trajectory where desorption is occurring.

Use `surface_normal` to define the direction "away" from the surface. Most commonly, this would be in positive z direction.

Use `extra_frames` to save additional steps before the desorption event begins.

A desorption is detected if the centre of mass of the molecule defined with `indices` is above `surface_distance_threshold` from the closest surface atom.
This is calculated with respect to `surface_normal` and will take into account periodic boundary conditions.
"""
OutputDesorptionTrajectory(
    indices;
    surface_normal = [0, 0, 1],
    surface_distance_threshold = austrip(4.0u"Å"),
    extra_frames = 0,
) = OutputDesorptionTrajectory(
    indices,
    convert(Vector{Float64}, surface_normal),
    surface_distance_threshold,
    extra_frames,
)
"""
    (output::OutputDesorptionTrajectory)(sol, i)

Only output parts of the trajectory where desorption is occurring.
"""
function (output::OutputDesorptionTrajectory)(sol, i)
    desorption_frame = Analysis.Diatomic.get_desorption_frame(
        sol.u,
        output.indices,
        sol.prob.p;
        surface_distance_threshold = output.surface_distance_threshold,
        surface_normal = output.surface_normal,
    )
    if isnothing(desorption_frame)
        return nothing
    end
    start_save_frame = desorption_frame - output.extra_frames
    if start_save_frame < 1
        return sol.u
    else
        return sol.u[start_save_frame:desorption_frame]
    end
end
export OutputDesorptionTrajectory

struct OutputDesorptionSnapshot{I<:Vector{Int},N<:Vector{Float64},D}
    indices::I
    surface_normal::N
    surface_distance_threshold::D
end
"""
    OutputDesorptionSnapshot(indices; surface_normal=[0, 0, 1], surface_distance_threshold=austrip(5.0u"Å"))

Save DynamicsVariables where desorption starts. (Same conditions as for OutputDesorptionTrajectory)

Use `surface_normal` to define the direction "away" from the surface. Most commonly, this would be in positive z direction.

Use `extra_frames` to save additional steps before the desorption event begins.

A desorption is detected if the centre of mass of the molecule defined with `indices` is above `surface_distance_threshold` from the closest surface atom.
This is calculated with respect to `surface_normal` and will take into account periodic boundary conditions.
"""
OutputDesorptionSnapshot(
    indices;
    surface_normal = [0, 0, 1],
    surface_distance_threshold = austrip(5.0u"Å"),
) = OutputDesorptionSnapshot(
    indices,
    convert(Vector{Float64}, surface_normal),
    surface_distance_threshold,
)
"""
    (output::OutputDesorptionTrajectory)(sol, i)

Only output trajectory snapshot where desorption begins. (Centre of mass velocity projected onto surface
normal changes sign)
"""
function (output::OutputDesorptionSnapshot)(sol, i)
    desorption_frame = Analysis.Diatomic.get_desorption_frame(
        sol.u,
        output.indices,
        sol.prob.p;
        surface_distance_threshold = output.surface_distance_threshold,
        surface_normal = output.surface_normal,
    )
    return isnothing(desorption_frame) ? nothing : sol.u[desorption_frame]
end
export OutputDesorptionSnapshot

"""
Outputs the instantaneous temperature of the selected atoms **in K** in the system at every save point.

Invoke with `OutputKineticTemperature(:)` for the entire system, or with `OutputKineticTemperature([1,2,3...])` for a subset of atoms.

"""
struct OutputKineticTemperature
    indices::Any
    OutputKineticTemperature(indices) = new(indices)
end

function (out::OutputKineticTemperature)(sol, i)
    # Allocate output vector
    kinetic_energies = zeros(
        typeof(DynamicsUtils.classical_kinetic_energy(sol.prob.p, sol.u[1])),
        length(sol.u),
    )
    # Determine number of atoms
    if isa(out.indices, Colon)
        n_atoms = length(sol.prob.p.atoms.masses)
    elseif isa(out.indices, UnitRange)
        n_atoms = length(collect(out.indices))
    else
        n_atoms = length(out.indices)
    end
    # Calculate kinetic temperatures.
    for snapshot in eachindex(sol.u)
        kinetic_energy = DynamicsUtils.classical_kinetic_energy(
            sol.prob.p.atoms.masses[out.indices],
            DynamicsUtils.get_velocities(sol.u[snapshot])[:, out.indices],
        )
        kinetic_energies[snapshot] = ustrip(
            uconvert(u"K", 2 * kinetic_energy * u"hartree/k_au") / ndofs(sol.prob.p) /
            n_atoms,
        )
    end
    return kinetic_energies
end

export OutputKineticTemperature

"""
    OutputNoise(sol, i)

Outputs the noise generated by the integrator at each time step in the trajectory.

**Note:** This requires using the `run_dynamics` command with `save_noise=true`, otherwise noise won't be accessible to this output function.
"""
function OutputNoise(sol, i)
    return sol.W
end

export OutputNoise

"""
    OutputEverything(sol,i)

Outputs the full DifferentialEquations solution object.

Storing this to disk is inefficient, but allows for full post-processing with
any of the functions defined in this module.
"""
OutputEverything(sol, i) = sol
OutputSolution(sol, i) = sol
export OutputEverything
export OutputSolution

end # module
