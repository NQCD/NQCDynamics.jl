
using AdvancedMH: DensityModel, RWMH, sample
using Distributions: MvNormal

function get_density_function(sim::Simulation, shape, temperature)
    density(R) = -evaluate_potential_energy(sim, reshape(R, shape)) / temperature
end

function get_density_function(sim::RingPolymerSimulation, shape, temperature)
    function density(R)
        transform_from_normal_modes!(sim.beads, R)
        -evaluate_potential_energy(sim, reshape(R, shape)) / temperature
    end
end

function sample_configurations(
    sim::AbstractSimulation,
    R0::AbstractArray,
    steps::Real;
    σ::AbstractArray=fill(0.1, size(R0)),
    )

    shape = size(R0)
    temperature = get_temperature(sim)
    density = get_density_function(sim, shape, temperature)

    density_model = DensityModel(density)
    d = MvNormal(σ[:])
    sampler = RWMH(d)

    chain = sample(density_model, sampler, convert(Int, steps); init_params=R0[:])

    [reshape(p.params, shape) for p in chain]
end
