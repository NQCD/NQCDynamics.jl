using NonadiabaticMolecularDynamics
using DifferentialEquations
using Plots
using Random

atoms = Atoms{Float64}([:H, :C])
cell = InfiniteCell{Float64}()
model = Models.Free()
sim = RingPolymerSimulation(1, 1e-2, cell, atoms, model, Dynamics.Classical(), 3, [:H])

R = cat([randn(1, 2) for i=1:3]..., dims=3)
P = cat([zeros(1, 2) for i=1:3]..., dims=3)
z = RingPolymerPhasespace(R, P)

problem = ODEProblem(Dynamics.motion!, z, (0.0, 1e3), sim)
solution = solve(problem)
plot(solution, vars=1:6)
