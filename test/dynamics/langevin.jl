push!(LOAD_PATH, pwd())
using Test
using NonadiabaticMolecularDynamics
using Unitful

@test Langevin{Float64}(1.0) == Langevin(1)

atomic_parameters = Atoms.AtomicParameters(Atoms.InfiniteCell(), [:H])
System{Langevin}(atomic_parameters, Models.Analytic.Free(), 10u"K", 1)

