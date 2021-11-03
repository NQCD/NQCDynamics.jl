using Test
using SafeTestsets
using NonadiabaticMolecularDynamics

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Core"
    @time @safetestset "Calculator Tests" begin include("Core/calculators.jl") end
    @time @safetestset "Simulation Tests" begin include("Core/simulations.jl") end
    @time @safetestset "Ring Polymer Tests" begin include("Core/ring_polymers.jl") end
    @time @safetestset "RingPolymerArrays Tests" begin include("Core/ring_polymer_array.jl") end
    @time @safetestset "Estimator tests" begin include("Core/estimators.jl") end
end

if GROUP == "All" || GROUP == "InitialConditions"
    @time @safetestset "Monte Carlo Tests" begin include("InitialConditions/monte_carlo.jl") end
    @time @safetestset "AdvancedMH Sampling Tests" begin include("InitialConditions/advancedmh_sampling.jl") end
    @time @safetestset "QuantisedDiatomic Tests" begin include("InitialConditions/quantised_diatomic.jl") end
    @time @safetestset "NonadiabaticDistributions" begin include("NonadiabaticDistributions/nonadiabatic_distributions.jl") end
    @time @safetestset "Distribution Tests" begin include("NonadiabaticDistributions/nuclear_distributions.jl") end
    @time @safetestset "Harmonic Wigner distribution tests" begin include("NonadiabaticDistributions/harmonic_wigner.jl") end
    @time @safetestset "ElectronicDistributions" begin include("NonadiabaticDistributions/electronic_distributions.jl") end
end

if GROUP == "All" || GROUP == "Dynamics"
    @time @safetestset "Classical Tests" begin include("Dynamics/classical.jl") end
    @time @safetestset "Langevin Tests" begin include("Dynamics/langevin.jl") end
    @time @safetestset "MDEF BAOAB Tests" begin include("Dynamics/mdef_baoab.jl") end
    @time @safetestset "MDEF Tests" begin include("Dynamics/mdef.jl") end
    @time @safetestset "RPMDEF Tests" begin include("Dynamics/rpmdef.jl") end
    @time @safetestset "BCBwithTsit5 Tests" begin include("Dynamics/bcbwithtsit5.jl") end
    @time @safetestset "FSSH Tests" begin include("Dynamics/fssh.jl") end
    @time @safetestset "NRPMD Tests" begin include("Dynamics/nrpmd.jl") end
    @time @safetestset "CMM Tests" begin include("Dynamics/cmm.jl") end
    @time @safetestset "Cell Boundary Callback Tests" begin include("Dynamics/cell_boundary_callback.jl") end
    @time @safetestset "Ehrenfest Tests" begin include("Dynamics/ehrenfest.jl") end
    @time @safetestset "IESH Tests" begin include("Dynamics/iesh.jl") end
    @time @safetestset "Ensemble Tests" begin include("Ensembles/ensembles.jl") end
    @time @safetestset "Ensemble Output Tests" begin include("Ensembles/outputs.jl") end
end
