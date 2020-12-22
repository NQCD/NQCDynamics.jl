using NonadiabaticMolecularDynamics
using DifferentialEquations
using Plots
using Unitful

atoms = Atoms{Float64}([:H, :H])
cell = InfiniteCell{Float64}()
calc = Calculators.Harmonic()
sim = Simulation(1, 300u"K", cell, atoms, calc, Dynamics.Langevin(1))

z = Dynamics.Phasespace(zeros(1, length(atoms)), zeros(1, length(atoms)))

problem = SDEProblem(Dynamics.motion!, Dynamics.random_force!, z, (0.0, 50), sim)
solution = solve(problem, EM(), dt=1e-2)
plot(solution, vars=range(atoms))