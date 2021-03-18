using SafeTestsets
using Test

@safetestset "Classical Tests" begin include("classical.jl") end
@safetestset "Langevin Tests" begin include("langevin.jl") end
@safetestset "MDEF BAOAB Tests" begin include("mdef_baoab.jl") end
@safetestset "MDEF Tests" begin include("mdef.jl") end
@safetestset "FSSH Tests" begin include("fssh.jl") end
# @safetestset "Fermionic bath ring polymer (experimental) Tests" begin include("fermionic.jl") end
@safetestset "NRPMD Tests" begin include("nrpmd.jl") end
@safetestset "Ensembles Tests" begin include("ensembles.jl") end
@safetestset "Cell Boundary Callback Tests" begin include("cell_boundary_callback.jl") end
