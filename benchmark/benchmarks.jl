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
    sim = Simulation{eval(type)}(atoms, model)
    z = DynamicsVariables(sim, v, r, SingleState(1))
    SUITE["TullyModelOne"]["classical"][type] = @benchmarkable run_trajectory($z, $tspan, $sim; dt=1.0)
end

types = [:NRPMD, :FSSH, :Ehrenfest]
n_beads = 10
r = fill(-5.0, 1, 1, n_beads)
v = fill(8.9, 1, 1, n_beads) ./ atoms.masses[1]
for type in types
    sim = RingPolymerSimulation{eval(type)}(atoms, model, n_beads)
    z = DynamicsVariables(sim, v, r, SingleState(1))
    SUITE["TullyModelOne"]["ring polymer"][type] = @benchmarkable run_trajectory($z, $tspan, $sim; dt=1.0)
end
