push!(LOAD_PATH, pwd())

using Test
using MDPrePostProcessing.IO

@test_nowarn read_system("test/io/slab.xyz")
