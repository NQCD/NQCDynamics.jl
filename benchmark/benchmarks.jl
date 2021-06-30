using NonadiabaticMolecularDynamics
using BenchmarkTools

const SUITE = BenchmarkGroup()
SUITE["TullyModelOne"] = BenchmarkGroup()
SUITE["TullyModelOne"]["classical"] = BenchmarkGroup()
SUITE["TullyModelOne"]["ring polymer"] = BenchmarkGroup()

model = TullyModelOne()
atoms = Atoms(2000)
types = [:FSSH, :Ehrenfest]

tspan = (0.0, 2500.0)

r = hcat(-5.0)
v = hcat(8.9) ./ atoms.masses[1]
for type in types
    sim = Simulation{eval(type)}(atoms, model; DoFs=1)
    z = Ensembles.select_u0(sim, v, r, 1, :diabatic)
    SUITE["TullyModelOne"]["classical"][type] = @benchmarkable Dynamics.run_trajectory($z, $tspan, $sim)
end

types = [:NRPMD, :FSSH, :Ehrenfest]
n_beads = 10
r = fill(-5.0, 1, 1, n_beads)
v = fill(8.9, 1, 1, n_beads) ./ atoms.masses[1]
for type in types
    sim = RingPolymerSimulation{eval(type)}(atoms, model, n_beads; DoFs=1)
    z = Ensembles.select_u0(sim, v, r, 1, :diabatic)
    SUITE["TullyModelOne"]["ring polymer"][type] = @benchmarkable Dynamics.run_trajectory($z, $tspan, $sim)
end

SUITE["DebyeSpinBoson"] = BenchmarkGroup()
SUITE["DebyeSpinBoson"]["classical"] = BenchmarkGroup()
SUITE["DebyeSpinBoson"]["ring polymer"] = BenchmarkGroup()

model = DebyeSpinBoson(50)
atoms = Atoms(ones(50))
types = [:FSSH, :Ehrenfest]
tspan = (0.0, 10.0)
r = randn(1, 50)
v = zeros(1, 50)
for type in types
    sim = Simulation{eval(type)}(atoms, model; DoFs=1)
    z = Ensembles.select_u0(sim, v, r, 1, :diabatic)
    SUITE["DebyeSpinBoson"]["classical"][type] = @benchmarkable Dynamics.run_trajectory($z, $tspan, $sim)
end

types = [:NRPMD, :FSSH, :Ehrenfest]
n_beads = 10
r = zeros(1, 50, n_beads)
v = randn(1, 50, n_beads)
for type in types
    sim = RingPolymerSimulation{eval(type)}(atoms, model, n_beads; DoFs=1)
    z = Ensembles.select_u0(sim, v, r, 1, :diabatic)
    SUITE["DebyeSpinBoson"]["ring polymer"][type] = @benchmarkable Dynamics.run_trajectory($z, $tspan, $sim)
end
