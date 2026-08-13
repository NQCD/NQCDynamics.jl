@testset "ASE" begin
    include("test_utils.jl")
    using PythonCall
    using NQCDInterfASE

    ase = pyimport("ase")

    h2 = ase.Atoms("H2", [0 0 0; 0 0 0.74])
    h2.center(vacuum=2.5)

    @testset "EMT" begin
        emt = pyimport("ase.calculators.emt")
        h2.calc = emt.EMT()
        model = AdiabaticASEModel(h2)
        @test test_model(model, 2)
    end

    @testset "GPAW" begin
        gpaw = pyimport("gpaw")
        h2.calc = gpaw.GPAW(xc="PBE", mode=gpaw.PW(300), txt="h2.txt")
        model = AdiabaticASEModel(h2)
        @test test_model(model, 2, rtol=1e-3)
    end

end
