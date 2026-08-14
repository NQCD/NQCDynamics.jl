using NQCDynamics, CUDA, NQCModels
using CairoMakie

struct ModelButGPU{M} <: NQCModels.ClassicalModels.ClassicalModel
    model::M
end

NQCModels.potential(model::ModelButGPU{M}, x) where M = NQCModels.potential(M, x) |> cu
function NQCModels.derivative!(model::ModelButGPU{M}, D, x) where M
    NQCModels.derivative!(M, D, x |> cu)
end


sim_langevin = Simulation{MDEF}(
    Atoms([20.0]),
    CompositeFrictionModel(
        Harmonic(m = 20.0),
        ConstantFriction(
            1,
            rand(1,1),
        ),
    ),
    temperature= 0,
)
u0_langevin_cpu = DynamicsVariables(sim_langevin, hcat(1.0), hcat(1.0))

@info "Running MDEF with 1D Harmonic + ConstantFriction on CPU"
langevin_cpu = run_dynamics(
    sim_langevin,
    (0, 2000.0),
    u0_langevin_cpu,
    trajectories= 1,
    output = (OutputPosition, OutputVelocity),
)

sim_langevin = Simulation{MDEF}(
    Atoms([20.0], NQCBase.CUDADevice),
    CompositeFrictionModel(
        Harmonic(m = 20.0),
        ConstantFriction(
            1,
            rand(1,1) |> cu,
        ),
    ),
    temperature= 0,
)
u0_langevin_cuda = DynamicsVariables(sim_langevin, hcat(1.0) |> cu, hcat(1.0) |> cu)

@info "Running MDEF with 1D Harmonic + ConstantFriction on GPU"
langevin_cuda = run_dynamics(
    sim_langevin,
    (0, 2000.0),
    u0_langevin_cuda,
    trajectories= 1,
    output = (OutputPosition, OutputVelocity),
)
