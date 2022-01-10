"""
    QuantisedDiatomic

This module exports two user facing functions:

- `generate_configurations`
    Creates a set of velocities and positions for diatomic molecule with specified
    vibrational `ν` and rotational `J` quantum numbers.

- `quantise_diatomic`
    Obtains vibrational `ν` and rotational `J` quantum numbers for a diatomic molecule
    with a given set of velocities and positions.

The central concept of this module is the EBK procedure which is nicely detailed here:
[Larkoski2006](@cite)

Inspired by VENUS96: [Hase1996](@cite)
"""
module QuantisedDiatomic

using QuadGK: QuadGK
using Rotations: Rotations
using LsqFit: LsqFit
using LinearAlgebra: norm, normalize!, nullspace, cross, I
using UnitfulAtomic: austrip
using Distributions: Uniform
using UnicodePlots: lineplot, lineplot!, DotCanvas
using Optim: Optim

using NQCDynamics: Simulation, Calculators, DynamicsUtils
using NQCBase: Atoms, PeriodicCell, InfiniteCell
using NQCModels: NQCModels
using NQCModels.AdiabaticModels: AdiabaticModel

export generate_configurations
export quantise_diatomic

struct SurfaceParameters{T}
    reduced_mass::T
    atom_indices::Vector{Int}
    slab::Matrix{T}
    height::T
    surface_normal::Vector{T}
    orthogonal_vectors::Matrix{T}
end

function SurfaceParameters(masses::AbstractVector, atom_indices::AbstractVector, slab::Matrix, height, surface_normal)
    μ = reduced_mass([masses[i] for i in atom_indices])
    orthogonal_vectors = nullspace(reshape(surface_normal, 1, :)) # Plane parallel to surface
    if size(slab, 2) == 0
        surface_top = 0
    else
        surface_top = maximum(slab[3,:])
    end
    return SurfaceParameters(μ, atom_indices, slab, height + surface_top,
        surface_normal, orthogonal_vectors)
end

struct GenerationParameters{T}
    direction::Vector{T}
    translational_energy::T
    samples::Int
end

Base.broadcastable(p::SurfaceParameters) = Ref(p)
Base.broadcastable(p::GenerationParameters) = Ref(p)

"""
    generate_configurations(sim, ν, J; samples=1000, height=10, normal_vector=[0, 0, 1],
        translational_energy=0, direction=[0, 0, -1], position=[0, 0, height])

Generate positions and momenta for given quantum numbers

`translational_energy`, `direction` and `position` specify the kinetic energy in
a specific direction with the molecule placed with centre of mass at `position`.

Keyword arguments `height` and `normal_vector` become relevant if the potential
requires specific placement of the molecule.
These allow the molecule to be placed at a distance `height` in the direction
`normal_vector` when performing potential evaluations.
"""
function generate_configurations(sim, ν, J;
    samples=1000,
    height=10.0,
    surface_normal=[0, 0, 1.0],
    translational_energy=0.0,
    direction=[0, 0, -1.0],
    atom_indices=[1, 2],
    r=zeros(size(sim)))

    _, slab = separate_slab_and_molecule(atom_indices, r)
    surface = SurfaceParameters(sim.atoms.masses, atom_indices, slab, height, surface_normal)
    generation = GenerationParameters(direction, austrip(translational_energy), samples)

    total_energy, bounds = find_total_energy(sim.calculator.model, ν, J, surface)

    radial_momentum = get_radial_momentum_function(total_energy, surface, J, sim.calculator.model)
    bonds, momenta = select_random_bond_lengths(bounds, radial_momentum, samples)

    @info "Generating the requested configurations..."
    configure_diatomic.(sim, bonds, momenta, J, surface, generation)
end

"""
    separate_slab_and_molecule(atom_indices, r)

Get the coordinates of the molecule and slab separately.
"""
function separate_slab_and_molecule(atom_indices, r)
    molecule = r[:,atom_indices]
    slab_indices = [i for i in axes(r, 2) if i ∉ atom_indices]
    slab = r[:,slab_indices]
    return (molecule, slab)
end

"""
    combine_slab_and_molecule(atom_indices, molecule, slab)

Revert the transformation `separate_slab_and_molecule`
"""
function combine_slab_and_molecule(atom_indices, molecule, slab)
    r = zeros(3, size(molecule, 2) + size(slab, 2))
    slab_indices = [i for i in axes(r, 2) if i ∉ atom_indices]
    r[:,atom_indices] .= molecule
    r[:,slab_indices] .= slab
    return r
end

function combine_slab_and_molecule(surface::SurfaceParameters, molecule)
    combine_slab_and_molecule(surface.atom_indices, molecule, surface.slab)
end

"""
    find_total_energy(sim, ν, J, surface)

Returns the energy associated with the specified quantum numbers
"""
function find_total_energy(model::AdiabaticModel, ν, J, surface::SurfaceParameters)

    μ = surface.reduced_mass

    bond_lengths = 0.5:0.01:5.0
    potential = calculate_diatomic_energy.(model, bond_lengths, surface) # Calculate binding curve
    potential_minimum, index = findmin(potential)
    r₀ = bond_lengths[index]
    potential .-= potential_minimum # Shift minimum to zero.

    k = calculate_force_constant(bond_lengths, potential)
    E = guess_initial_energy(k, μ, ν, potential_minimum)

    V(r) = effective_potential(r, J, model, surface)

    "Evaluate vibrational quantum number."
    function nᵣ(E, r₁, r₂)
        kernel(r) = sqrt(E - V(r))
        integral, _ = QuadGK.quadgk(kernel, r₁, r₂; maxevals=100)
        return sqrt(2μ)/π * integral - 1/2
    end

    "Derivative of the above"
    function ∂nᵣ∂E(E, r₁, r₂)
        kernel(r) = 1 / sqrt(E - V(r))
        integral, _ = QuadGK.quadgk(kernel, r₁, r₂; maxevals=100)
        return sqrt(2μ)/2π * integral
    end

    "Second derivative of the above"
    function ∂²nᵣ∂E²(E, r₁, r₂)
        kernel(r) = 1 / sqrt(E - V(r))^3
        integral, _ = QuadGK.quadgk(kernel, r₁, r₂; maxevals=100)
        return -sqrt(2μ)/4π * integral
    end

    """
    Evaluates the function, gradient and hessian as described in the Optim documentation.
    """
    function fgh!(F,G,H,x)
        E = x[1]
        r₁, r₂ = find_integral_bounds(model, E, J, surface, r₀)

        νtmp = nᵣ(E, r₁, r₂)
        difference = νtmp - ν

        if G !== nothing
            G[1] = ∂nᵣ∂E(E, r₁, r₂) * sign(difference)
        end
        if H !== nothing
            H[1] = ∂²nᵣ∂E²(E, r₁, r₂) * sign(difference)
        end
        if F !== nothing
            return abs(difference)
        end
    end

    @info "Converging the total energy to match the chosen quantum numbers..."
    optim = Optim.optimize(Optim.only_fgh!(fgh!), [E], Optim.Newton(),
        Optim.Options(show_trace=true, x_tol=1e-3)
    )
    E = Optim.minimizer(optim)[1]
    bounds = find_integral_bounds(model, E, J, surface, r₀)

    E, bounds
end

"""
    quantise_diatomic(sim::Simulation, v::Matrix, r::Matrix;
        height=10, normal_vector=[0, 0, 1])

Quantise the vibrational and rotational degrees of freedom for the specified
positions and velocities

When evaluating the potential, the molecule is moved to `height` in direction `normal_vector`.
If the potential is independent of centre of mass position, this has no effect.
Otherwise, be sure to modify these parameters to give the intended behaviour.
"""
function quantise_diatomic(sim::Simulation, v::Matrix, r::Matrix;
    height=10.0, surface_normal=[0, 0, 1.0], atom_indices=[1, 2])


    r, slab = separate_slab_and_molecule(atom_indices, r)
    surface = SurfaceParameters(sim.atoms.masses, atom_indices, slab, austrip(height), surface_normal)

    r_com = subtract_centre_of_mass(r, [sim.atoms.masses[i] for i in atom_indices])
    v_com = subtract_centre_of_mass(v, [sim.atoms.masses[i] for i in atom_indices])
    p_com = v_com .* [sim.atoms.masses[i] for i in atom_indices]'

    k = DynamicsUtils.classical_kinetic_energy(sim, v_com)
    p = calculate_diatomic_energy(sim.calculator.model, bond_length(r_com), surface)
    E = k + p

    L = total_angular_momentum(r_com, p_com)
    J = (sqrt(1+4*L^2) - 1) / 2 # L^2 = J(J+1)ħ^2

    r₁, r₂ = find_integral_bounds(sim.calculator.model, E, J, surface, bond_length(r_com))

    μ = surface.reduced_mass
    V(r) = effective_potential(r, J, sim.calculator.model, surface)

    function nᵣ(E, r₁, r₂)
        kernel(r) = sqrt(E - V(r))
        integral, _ = QuadGK.quadgk(kernel, r₁, r₂; maxevals=100)
        return sqrt(2μ)/π * integral - 1/2
    end

    ν = nᵣ(E, r₁, r₂)

    ν, J
end

rotational_energy(r, μ, J) = J*(J+1) / (2μ*r^2)

function effective_potential(r, J, model, surface)
    rotational = rotational_energy(r, surface.reduced_mass, J)
    potential = calculate_diatomic_energy(model, r, surface)
    return rotational + potential
end

function get_radial_momentum_function(total_energy, surface, J, model)
    radial_momentum(r) = sqrt(2surface.reduced_mass * (total_energy - effective_potential(r, J, model, surface)))
end

function find_integral_bounds(model, total_energy, J, surface, r₀, bond_limits=(0.5, 5.0))
    optim_func(r) = abs(total_energy - effective_potential(r, J, model, surface))
    optim = Optim.optimize(optim_func, first(bond_limits), r₀)
    r₁ = Optim.minimizer(optim)
    optim = Optim.optimize(optim_func, r₀, last(bond_limits))
    r₂ = Optim.minimizer(optim)
    return r₁, r₂
end

"""
Pick a random bond length and corresponding radial momentum that matches the
radial probability distribution.

Uses rejection sampling: https://en.wikipedia.org/wiki/Rejection_sampling
"""
function select_random_bond_lengths(bounds, momentum_function, samples)
    @info "Generating distribution of bond lengths and radial momenta..."
    M = 1 / momentum_function(bounds[1]+1e-3)
    radius_proposal = Uniform(bounds...)
    bonds = zeros(samples)
    momenta = zeros(samples)
    for i=1:samples
        keep_going = true
        while keep_going
            r = rand(radius_proposal) # Generate radius
            p = momentum_function(r) # Evaluate probability
            probability = 1 / p # Probability inversely proportional to momentum
            if rand() < probability / M
                bonds[i] = r
                momenta[i] = rand([-1, 1]) * p
                keep_going = false
            end
        end
    end
    bonds, momenta
end

"""
Randomly orient molecule in space for a given bond length and radial momentum
"""
function configure_diatomic(sim, bond, momentum, J, surface::SurfaceParameters, generation::GenerationParameters)
    r = zeros(3,2)
    v = zeros(3,2)
    r[1,1] = bond
    v[1,1] = momentum ./ surface.reduced_mass
    r = subtract_centre_of_mass(r, [sim.atoms.masses[i] for i in surface.atom_indices])
    v = subtract_centre_of_mass(v, [sim.atoms.masses[i] for i in surface.atom_indices])

    L = sqrt(J * (J + 1))

    random = 2π*rand()
    ω = L .* [0, sin(random), cos(random)] ./ (surface.reduced_mass * bond^2)

    v[:,1] .+= cross(ω, r[:,1])
    v[:,2] .+= cross(ω, r[:,2])

    apply_random_rotation!(v, r)
    position_above_surface!(r, surface.height, sim.cell)
    apply_translational_impulse!(v, [sim.atoms.masses[i] for i in surface.atom_indices], generation.translational_energy, generation.direction)

    v, r
end

"""
Randomly rotate each column of two 3*N matrix, same rotation for all columns.
"""
function apply_random_rotation!(x, y)
    rotation = rand(Rotations.RotMatrix{3})
    for (col1, col2) in zip(eachcol(x), eachcol(y))
        col1 .= rotation * col1
        col2 .= rotation * col2
    end
end

function position_above_surface!(r, height, cell::PeriodicCell)
    r[3,:] .+= height
    a1 = cell.vectors[:,1]
    a2 = cell.vectors[:,2]
    displacement = rand()*a1+rand()*a2
    r .+= displacement
    return nothing
end

function position_above_surface!(r, height, ::InfiniteCell)
    r[3,:] .+= height
    return nothing
end

function apply_translational_impulse!(v, masses, translational_energy, direction)
    velocity_impulse = velocity_from_energy(masses, translational_energy)
    v .+= velocity_impulse .* normalize!(direction)
end

function velocity_from_energy(masses, energy)
    m = sum(masses)
    sqrt(2austrip(energy) / m)
end

"""
    calculate_diatomic_energy(model::AdiabaticModel, bond_length::Real;
        height=10, normal_vector=[0, 0, 1])

Returns potential energy of diatomic with `bond_length` at `height` from surface.

Orients molecule parallel to the surface at the specified height, the surface is assumed
to intersect the origin.
This requires that the model implicitly provides the surface, or works fine without one.
"""
function calculate_diatomic_energy(
    model::AdiabaticModel, r₀::Real, surface::SurfaceParameters
)

    (;surface_normal, orthogonal_vectors, height) = surface
    molecule = surface_normal .* height .+ orthogonal_vectors .* r₀ ./ sqrt(2)
    R = combine_slab_and_molecule(surface, molecule)

    return NQCModels.potential(model, R)
end

"""
    calculate_force_constant(sim::Simulation, surface::SurfaceParameters)

Evaluate energy for different bond lengths and identify force constant. 

This uses LsqFit.jl to fit a morse potential with optimised minimum.
"""
function calculate_force_constant(bond_lengths, potential)

    morse!(f, x, p) = (@. f = p[1]*(1 - exp(-p[2]*(x-p[3])))^2)
    function jacobian!(J, x, p)
        @. J[:,1] = (1 - exp(-p[2]*(x-p[3])))^2
        @. J[:,2] = 2p[1]*(1 - exp(-p[2]*(x-p[3]))) * exp(-p[2]*(x-p[3])) * (x - p[3])
        @. J[:,3] = -2p[1]*(1 - exp(-p[2]*(x-p[3]))) * exp(-p[2]*(x-p[3])) * p[2]
    end

    σ = potential[end]
    weights = exp.(-potential.^2 ./ σ^2) # Weight the area around the minimum most heavily

    fit = LsqFit.curve_fit(morse!, jacobian!, bond_lengths, potential, weights, [1.0, 1.0, 1.0];
        lower=[0.0, 0.0, 0.0], inplace=true)

    plot_binding_curve(morse!, bond_lengths, potential, fit.param)

    force_constant = fit.param[2]^2 * 2 * fit.param[1]
    return force_constant
end

function plot_binding_curve(morse!, bond_lengths, V, parameters)
    fitted_points = zero(bond_lengths)
    morse!(fitted_points, bond_lengths, parameters)
    plt = lineplot(bond_lengths, V;
        title="Binding curve", xlabel="Bond length / bohr", ylabel="Energy / Hartree",
        name="Actual values", canvas=DotCanvas, border=:ascii
    )
    lineplot!(plt, bond_lengths, fitted_points; name="Fitted curve")
    show(plt)
    println()
    @info "The two lines shown above should closely match. \
        This indicates the evaluation of the potential is working correctly."
end

function guess_initial_energy(k, μ, ν, potential_minimum)
    ω = sqrt(k / μ)
    E_guess = (ν + 1/2) * ω
    return E_guess + potential_minimum
end

subtract_centre_of_mass(x, m) = x .- centre_of_mass(x, m)
@views centre_of_mass(x, m) = (x[:,1].*m[1] .+ x[:,2].*m[2]) ./ (m[1]+m[2])
reduced_mass(atoms::Atoms) = reduced_mass(atoms.masses)
reduced_mass(m::AbstractVector) = m[1]*m[2]/(m[1]+m[2])
@views bond_length(r) = norm(r[:,1] .- r[:,2])
@views total_angular_momentum(r, p) = norm(cross(r[:,1], p[:,1]) + cross(r[:,2], p[:,2]))

end # module
