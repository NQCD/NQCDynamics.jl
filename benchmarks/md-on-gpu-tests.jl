### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ dcff3f80-8513-11f1-8ee4-594e3df1267c
import Pkg

# ╔═╡ 17a1dcb4-7f1d-4e45-8892-21dd97e0c8fe
Pkg.activate(".")

# ╔═╡ 38d1ffd7-a731-4341-bb60-edcb2ecf0ce3
using NQCDynamics, CUDA, NQCModels

# ╔═╡ 7c267a8d-ff6b-480d-8dfb-4af47ba9e0d8
using CairoMakie

# ╔═╡ 68667991-e081-4702-92f4-d6db56e69df1
md"""
# Testing GPU-compatibility for different dynamics methods

| Dynamics Method | Model used | Integrator | Is working? |
|---|---|---|---|
| Classical | `DiatomicHarmonic` | `VelocityVerlet` | ✅ |
| Langevin | `Harmonic` | `StochasticDiffEq.BAOAB` | ❌ |

"""

# ╔═╡ b1f6de4e-7c73-43c9-a953-baa1f07c6dc4
md"# Classical dynamics on CUDA"

# ╔═╡ a59f4684-d877-4a46-be76-72fc464d4280


# ╔═╡ cc1d2fbc-f572-49c2-b98d-9277b358da29
begin # CPU dynamics example
    sim = Simulation(Atoms([1, 1]), DiatomicHarmonic())
    v = rand(3, 2)
    u0 = DynamicsVariables(sim, zeros(3, 2), hcat(randn(3), randn(3).+1))
    
    traj = run_dynamics(sim, (0.0, 10.0), u0; dt=0.1, output=OutputPosition, trajectories = 100)
end

# ╔═╡ 5abb828f-d525-41b0-aef0-8168ec4f88a8
begin # GPU dynamics
    u0_gpu = DynamicsVariables(sim, zeros(3, 2) |> cu, hcat(randn(3), randn(3).+1) |> cu)
    trj_gpu = run_dynamics(sim, (0.0, 10.0), u0; dt=0.1, output=OutputPosition, trajectories = 100)
end

# ╔═╡ e1be6d37-1a6c-4084-99d6-22ffb5bcd596
md"## Langevin dynamics"

# ╔═╡ 2e7fa650-2b47-4ff7-8bd6-030d0ea40766
md"""
- Constant friction. 
- Uses `StochasticDiffEq`'s `BAOAB` integrator. 
- Harmonic oscillator
"""

# ╔═╡ feb8c3d5-fc87-47a6-ade4-e82eea8dcc82
begin # Langevin example
    sim_langevin = Simulation{Langevin}(
        Atoms([20.0]),
        Harmonic(m = 20.0),
        γ = 0.02,
        temperature= 0,
    )
    u0_langevin_cpu = DynamicsVariables(sim_langevin, hcat(1.0), hcat(1.0))
end

# ╔═╡ 64ea24dc-c625-4047-826c-2b0e282c58ab
begin # Langevin GPU
    u0_langevin_cuda = DynamicsVariables(sim_langevin, CuArray(hcat(1.0)), CuArray(hcat(1.0)))
end

# ╔═╡ 9c565cf4-a95d-4063-855a-ba587f8cc90b
langevin_cpu = run_dynamics(
    sim_langevin,
    (0, 2000.0),
    u0_langevin_cpu,
    trajectories= 1,
    output = (OutputPosition, OutputVelocity),
)

# ╔═╡ 5ee502e7-6877-444b-be90-be411243f886
langevin_cuda = run_dynamics(
    sim_langevin,
    (0, 2000.0),
    u0_langevin_cuda,
    trajectories= 1,
    output = (OutputPosition, OutputVelocity),
)

# ╔═╡ e79d372e-fee6-4d22-9b8b-9e3c4191abbe
lines(langevin_cpu[:OutputPosition] .|> first, axis = (title = "Langevin Dynamics CPU", ), )

# ╔═╡ b53315bb-73ee-4dfc-adac-2dd9a66f9d6d
md"# Ehrenfest dynamics"

# ╔═╡ b0270e88-fdaa-4e28-ab51-43a19f6fd07e
ehrenfest_sim = Simulation{Ehrenfest}(
    Atoms(1980),
    AnanthModelOne(),
)

# ╔═╡ eef4d6d3-aaf0-4119-889f-daf8e4f61c37
ehrenfest_u0_cpu = DynamicalDistribution(
    sqrt(0.03 * 2 * ehrenfest_sim.atoms.masses[1]),
    NQCDynamics.Distributions.Normal(-5, 1/sqrt(0.25)),
    size(ehrenfest_sim),
) * PureState(1, Adiabatic())

# ╔═╡ d6940a24-eff1-40ba-98e0-cfb5232f4464
ehrenfest_dynamics_cpu = run_dynamics(
    ehrenfest_sim,
    (0.0,3000.0),
    ehrenfest_u0_cpu,
    trajectories= 10,
    output= (OutputVelocity, OutputPosition),
    dt = 1.0,
)

# ╔═╡ c16f1060-7619-4695-af9d-0e289d367a37
NQCDynamics.Distributions.Normal

# ╔═╡ 4a99f5e6-1ddc-4329-8b23-05a8642e095f
ehrenfest_u0_cuda = DynamicalDistribution(
    sqrt(0.03 * 2 * ehrenfest_sim.atoms.masses[1]) |> hcat |> CuArray,
    CUDA.rand(1,1),
    size(ehrenfest_sim),
) * PureState(1, Adiabatic())

# ╔═╡ 35af5fb7-19ed-49c6-bc0d-1687f049f93a
ehrenfest_dynamics_cuda = run_dynamics(
    ehrenfest_sim,
    (0.0,3000.0),
    ehrenfest_u0_cuda,
    trajectories= 10,
    output= (OutputVelocity, OutputPosition),
    dt = 1.0,
)

# ╔═╡ Cell order:
# ╠═dcff3f80-8513-11f1-8ee4-594e3df1267c
# ╠═17a1dcb4-7f1d-4e45-8892-21dd97e0c8fe
# ╠═38d1ffd7-a731-4341-bb60-edcb2ecf0ce3
# ╠═7c267a8d-ff6b-480d-8dfb-4af47ba9e0d8
# ╠═68667991-e081-4702-92f4-d6db56e69df1
# ╟─b1f6de4e-7c73-43c9-a953-baa1f07c6dc4
# ╟─a59f4684-d877-4a46-be76-72fc464d4280
# ╟─cc1d2fbc-f572-49c2-b98d-9277b358da29
# ╟─5abb828f-d525-41b0-aef0-8168ec4f88a8
# ╟─e1be6d37-1a6c-4084-99d6-22ffb5bcd596
# ╟─2e7fa650-2b47-4ff7-8bd6-030d0ea40766
# ╟─feb8c3d5-fc87-47a6-ade4-e82eea8dcc82
# ╟─64ea24dc-c625-4047-826c-2b0e282c58ab
# ╟─9c565cf4-a95d-4063-855a-ba587f8cc90b
# ╟─5ee502e7-6877-444b-be90-be411243f886
# ╠═e79d372e-fee6-4d22-9b8b-9e3c4191abbe
# ╠═b53315bb-73ee-4dfc-adac-2dd9a66f9d6d
# ╠═b0270e88-fdaa-4e28-ab51-43a19f6fd07e
# ╠═eef4d6d3-aaf0-4119-889f-daf8e4f61c37
# ╠═d6940a24-eff1-40ba-98e0-cfb5232f4464
# ╠═c16f1060-7619-4695-af9d-0e289d367a37
# ╠═4a99f5e6-1ddc-4329-8b23-05a8642e095f
# ╠═35af5fb7-19ed-49c6-bc0d-1687f049f93a
