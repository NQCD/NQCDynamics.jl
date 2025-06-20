using Test
using SafeTestsets
using NQCDynamics
using JSON

const GROUP = get(ENV, "GROUP", "All")

#= if GROUP == "All" || GROUP == "Core"
    @safetestset "FastDeterminant Tests" begin
        include("Core/FastDeterminant.jl")
    end
    @safetestset "Simulation Tests" begin
        include("Core/simulations.jl")
    end
    @safetestset "Ring Polymer Tests" begin
        include("Core/ring_polymers.jl")
    end

end

if GROUP == "All" || GROUP == "Analysis"
    @safetestset "Estimator tests" begin
        include("Core/estimators.jl")
    end
    @safetestset "Structure Manipulation" begin
        include("Structure/structure.jl")
    end
    @safetestset "Diatomic Analysis" begin
        include("Analysis/diatomic.jl")
    end
end

if GROUP == "All" || GROUP == "InitialConditions"
    @safetestset "Monte Carlo Tests" begin
        include("InitialConditions/monte_carlo.jl")
    end
    @safetestset "AdvancedMH Sampling Tests" begin
        include("InitialConditions/advancedmh_sampling.jl")
    end
    @safetestset "QuantisedDiatomic Tests" begin
        include("InitialConditions/quantised_diatomic.jl")
    end
end
 =#
if GROUP == "All" || GROUP == "dynamics_classical"
    @safetestset "DynamicsUtils Tests" begin
        include("Dynamics/DynamicsUtils.jl")
    end
    @safetestset "Classical Tests" begin
        include("Dynamics/classical.jl")
    end
    @safetestset "BCB Tests" begin
        include("Dynamics/algorithms/bcb.jl")
    end
    @safetestset "Langevin Tests" begin
        include("Dynamics/langevin.jl")
    end
end

if GROUP == "All" || GROUP == "dynamics_mdef"
    @safetestset "MDEF BAOAB Tests" begin
        include("Dynamics/mdef_baoab.jl")
    end
    @safetestset "MDEF Tests" begin
        include("Dynamics/mdef.jl")
    end
    @safetestset "DiabaticMDEF Tests" begin
        include("Dynamics/diabatic_mdef.jl")
    end
    @safetestset "RPMDEF Tests" begin
        include("Dynamics/rpmdef.jl")
    end
end

if GROUP == "All" || GROUP == "dynamics_surface_hopping"
    @safetestset "BCBwithTsit5 Tests" begin
        include("Dynamics/bcbwithtsit5.jl")
    end
    @safetestset "FSSH Tests" begin
        include("Dynamics/fssh.jl")
    end
     @safetestset "IESH Tests" begin
        include("Dynamics/iesh.jl")
    end
    @safetestset "Decoherence Tests" begin
        include("Dynamics/test_decoherence_corrections.jl")
    end
    @safetestset "RPIESH Tests" begin
        include("Dynamics/rpiesh.jl")
    end
end

if GROUP == "All" || GROUP == "dynamics_cme"
    @safetestset "CME Tests" begin
        include("Dynamics/cme.jl")
    end
    @safetestset "BCME Tests" begin
        include("Dynamics/bcme.jl")
    end
    @safetestset "RPCME Tests" begin
        include("Dynamics/rpcme.jl")
    end
    @safetestset "RPBCME Tests" begin
        include("Dynamics/rpbcme.jl")
    end
end


if GROUP == "All" || GROUP == "dynamics_mapping"
    @safetestset "Electronic dynamics Tests" begin
        include("Dynamics/electronic_dynamics.jl")
    end
    @safetestset "NRPMD Tests" begin
        include("Dynamics/nrpmd.jl")
    end
    @safetestset "CMM Tests" begin
        include("Dynamics/cmm.jl")
    end
    @safetestset "RPeCMM Tests" begin
        include("Dynamics/rpecmm.jl")
    end
    @safetestset "SpinMapping Tests" begin
        include("Dynamics/spin_mapping.jl")
    end
end

if GROUP == "All" || GROUP == "dynamics_ehrenfest"
    @safetestset "Ehrenfest Tests" begin
        include("Dynamics/ehrenfest.jl")
    end
    @safetestset "NA Ehrenfest Tests" begin
        include("Dynamics/ehrenfest_na.jl")
    end
    @safetestset "NA RPEhrenfest Tests" begin
        include("Dynamics/rp_ehrenfest_na.jl")
    end
end


if GROUP == "All" || GROUP == "dynamics_ensembles"
    @safetestset "Ensemble Tests" begin
        include("Ensembles/ensembles.jl")
    end
    @safetestset "Ensemble Selection Tests" begin
        include("Ensembles/selections.jl")
    end
    @safetestset "Ensemble Reduction Tests" begin
        include("Ensembles/reductions.jl")
    end
    @safetestset "Ensemble Output Tests" begin
        include("Ensembles/outputs.jl")
    end
    @safetestset "Cell Boundary Callback Tests" begin
        include("Dynamics/cell_boundary_callback.jl")
    end
end

