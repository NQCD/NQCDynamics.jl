using Test
using NQCModels

m = ErpenbeckThoss(;Γ=0.1)
b = TrapezoidalRule(10, -10, 10)
@testset "nelectrons" begin
    model = AndersonHolstein(m, b) # fermilevel defaults to zero
    @test NQCModels.nelectrons(model) == 5

    for (n, f) in zip([1, 6, 10], [-10, 2, 10])
        model = AndersonHolstein(m, b; fermi_level=f)
        @test NQCModels.nelectrons(model) == n
    end
end
