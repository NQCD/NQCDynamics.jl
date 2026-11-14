
struct ProductDistribution{V,R,E<:ElectronicDistribution}
    nuclear::DynamicalDistribution{V,R}
    electronic::E
end

Base.:*(N::DynamicalDistribution, E::ElectronicDistribution) = ProductDistribution(N, E)
Base.:*(E::ElectronicDistribution, N::DynamicalDistribution) = ProductDistribution(N, E)
