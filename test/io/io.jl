using Test
using NonadiabaticMolecularDynamics.IO

@test_nowarn read_system("test/io/slab.xyz")
