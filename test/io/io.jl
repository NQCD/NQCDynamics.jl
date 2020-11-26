push!(LOAD_PATH, pwd())

using Test
using MDPrePostProcessing.Basic

@test_nowarn p, z = read_system("test/io/slab.xyz")
