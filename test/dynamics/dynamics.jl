using SafeTestsets
using Test

@time @safetestset "Classical Tests" begin include("classical.jl") end
@time @safetestset "Langevin Tests" begin include("langevin.jl") end
@time @safetestset "MDEF BAOAB Tests" begin include("mdef_baoab.jl") end
@time @safetestset "MDEF Tests" begin include("mdef.jl") end
@time @safetestset "FSSH Tests" begin include("fssh.jl") end
# @safetestset "Fermionic bath ring polymer (experimental) Tests" begin include("fermionic.jl") end
@time @safetestset "NRPMD Tests" begin include("nrpmd.jl") end
@time @safetestset "Cell Boundary Callback Tests" begin include("cell_boundary_callback.jl") end
