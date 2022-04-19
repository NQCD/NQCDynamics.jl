using AdvancedMH: MetropolisHastings, DensityModel, sample, StaticProposal

struct GenerationParameters{T}
    direction::Vector{T}
    translational_energy::T
    samples::Int
end

Base.broadcastable(p::GenerationParameters) = Ref(p)

"""
Pick a random bond length and corresponding radial momentum that matches the
radial probability distribution.
"""
function select_random_bond_lengths(bounds, momentum_function, samples)
    @info "Generating distribution of bond lengths and radial momenta..."

    insupport(θ) = bounds[1] < θ < bounds[2]
    density(θ) = insupport(θ) ? log(1/momentum_function(θ)) : -Inf
    model = DensityModel(density)
    spl = MetropolisHastings(StaticProposal(Uniform(bounds...)))
    chain = sample(model, spl, samples; chain_type=Vector{NamedTuple})

    bonds = [c.param_1 for c in chain]
    momenta = [1/exp(c.lp) * rand((-1, 1)) for c in chain]

    return bonds, momenta
end

function plot_distributions(bonds, velocities)
    r_hist = histogram(bonds, title="Bond length distribution", border=:ascii)
    v_hist = histogram(velocities, title="Radial velocity distribution", border=:ascii)
    show(r_hist)
    println()
    show(v_hist)
    println()
end

"""
Randomly orient molecule in space for a given bond length and radial momentum
"""
function configure_diatomic(sim, bond, velocity, J, environment::EvaluationEnvironment, generation::GenerationParameters, μ)
    r = zeros(3,2)
    v = zeros(3,2)
    r[1,1] = bond
    v[1,1] = velocity
    r = subtract_centre_of_mass(r, masses(sim)[environment.molecule_indices])
    v = subtract_centre_of_mass(v, masses(sim)[environment.molecule_indices])

    L = sqrt(J * (J + 1))

    random = 2π*rand()
    ω = L .* [0, sin(random), cos(random)] ./ (μ * bond^2)

    v[:,1] .+= cross(ω, r[:,1])
    v[:,2] .+= cross(ω, r[:,2])

    apply_random_rotation!(v, r)
    position_above_surface!(r, environment.offset, sim.cell)
    apply_translational_impulse!(v, masses(sim)[environment.molecule_indices], generation.translational_energy, generation.direction)

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

function position_above_surface!(r::AbstractMatrix, offset::AbstractVector, cell::PeriodicCell)
    r .+= offset
    a1 = cell.vectors[:,1]
    a2 = cell.vectors[:,2]
    displacement = rand()*a1+rand()*a2
    r .+= displacement
    return nothing
end

function position_above_surface!(r::AbstractMatrix, offset::AbstractVector, ::InfiniteCell)
    r .+= offset
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