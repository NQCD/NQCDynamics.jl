using Test
using NQCDistributions
using Unitful, UnitfulAtomic
using DataInterpolations: LinearInterpolation

@testset "PureState" begin
    d = PureState(2)
    @test NQCDistributions.density_matrix(d, 2) == [0 0; 0 1]
end

@testset "MixedState" begin
    d = MixedState([0.3, 0.4])
    @test NQCDistributions.density_matrix(d) == [0.3 0; 0 0.4]
end

@testset "FermiDiracState" begin
    d = FermiDiracState(0.0, 300u"K")
    NQCDistributions.density_matrix(d, austrip.((-1:0.01:1) .* u"eV"))
end

@testset "density_matrix is exported" begin
    @test density_matrix(PureState(2), 2) == [0 0; 0 1]
end

@testset "NonEq-ElectronicDist" begin
    energies = collect(range(-1.0, 1.0; length=11))
    occupations = [ϵ < 0 ? 1.0 : 0.0 for ϵ in energies]
    dis_spline = LinearInterpolation(occupations, energies)
    dos_spline = LinearInterpolation(ones(length(energies)), energies)

    @testset "NonEqState" begin
        d = NonEqState(dis_spline, dos_spline)
        @test d isa NQCDistributions.ElectronicDistribution
        @test d.statetype == Adiabatic()
        @test d.available_states == Colon()
        @test d.dis_spline(-0.5) ≈ 1.0
        @test d.dis_spline(0.5) ≈ 0.0
        @test d.dos_spline(0.3) ≈ 1.0

        d = NonEqState(dis_spline, dos_spline; statetype=Diabatic(), available_states=1:5)
        @test d.statetype == Diabatic()
        @test d.available_states == 1:5
    end

    @testset "DensityMatrix" begin
        # Fermi-Dirac diagonal, Fermi level and temperature in atomic units
        d = DensityMatrix(0.0, 0.01)
        @test d.diagonal isa FermiDiracState
        @test d.diagonal.fermi_level == 0.0
        @test d.diagonal.β ≈ 1 / 0.01
        @test d.diagonal.statetype == Adiabatic()
        @test d.systemtype == Binary()

        d = DensityMatrix(0.0, 0.01; statetype=Diabatic(), systemtype=Statistical())
        @test d.diagonal.statetype == Diabatic()
        @test d.systemtype == Statistical()

        # Unlike `FermiDiracState`, the constructor only accepts `Float64`, so Unitful quantities fail.
        @test_broken DensityMatrix(0.0u"eV", 300u"K") isa DensityMatrix

        # Non-equilibrium diagonal
        d = DensityMatrix(dis_spline, dos_spline)
        @test d.diagonal isa NonEqState
        @test d.diagonal.dis_spline === dis_spline
        @test d.systemtype == Binary()
    end
end
