using NonadiabaticMolecularDynamics
using BenchmarkTools

const SUITE = BenchmarkGroup()
SUITE["dynamics"] = BenchmarkGroup()
SUITE["dynamics"]["classical"] = BenchmarkGroup()
SUITE["dynamics"]["ring polymer"] = BenchmarkGroup()

model = TullyModelOne()
atoms = Atoms(2000)
types = [:FSSH, :Ehrenfest]

tspan = (0.0, 2500.0)

r = hcat(-5.0)
v = hcat(8.9) ./ atoms.masses[1]
for type in types
    sim = Simulation{eval(type)}(atoms, model; DoFs=1)
    z = Ensembles.select_u0(sim, v, r, 1, :diabatic)
    SUITE["dynamics"]["classical"][type] = @benchmarkable Dynamics.run_trajectory($z, $tspan, $sim)
end

types = [:NRPMD, :FSSH, :Ehrenfest]

n_beads = 10
r = fill(-5.0, 1, 1, n_beads)
v = fill(8.9, 1, 1, n_beads) ./ atoms.masses[1]
for type in types
    sim = RingPolymerSimulation{eval(type)}(atoms, model, n_beads; DoFs=1)
    z = Ensembles.select_u0(sim, v, r, 1, :diabatic)
    SUITE["dynamics"]["ring polymer"][type] = @benchmarkable Dynamics.run_trajectory($z, $tspan, $sim)
end
