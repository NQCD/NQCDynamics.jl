using Test
using NQCDistributions
using Unitful, UnitfulAtomic

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
