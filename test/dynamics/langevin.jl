using Test
using NonadiabaticMolecularDynamics
using Unitful

@test Langevin{Float64}(1.0, 1.0) == Langevin(1, 1)
@test Langevin{Float64}(1.0, 1.0) == Langevin(1.0, 1)
@test_nowarn Langevin(1u"K", 1)

atomic_parameters = Atoms.AtomicParameters(InfiniteCell(), [:H])
System{Langevin}(atomic_parameters, Models.Analytic.Free(), 10u"K", 1)

