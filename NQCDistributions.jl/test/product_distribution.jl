using NQCDistributions

@test DynamicalDistribution(Normal(), Normal(), (3,2)) * PureState(10) isa NQCDistributions.ProductDistribution
