
struct GenerationParameters{T}
    direction::Vector{T}
    translational_energy::T
    samples::Int
end

Base.broadcastable(p::GenerationParameters) = Ref(p)

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
function configure_diatomic(sim, bond, momentum, J, environment::EvaluationEnvironment, generation::GenerationParameters, μ)
    r = zeros(3,2)
    v = zeros(3,2)
    r[1,1] = bond
    v[1,1] = momentum ./ μ
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