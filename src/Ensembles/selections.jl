
abstract type AbstractSelection end

"""
$(TYPEDEF)

Select the initial conditions from the distribution in order. 
"""
struct OrderedSelection{D} <: AbstractSelection
    "Distribution that is sampled."
    distribution::D
end

"""
$(TYPEDEF)

Obtain initial conditions by randomly sampling the distribution.
"""
struct RandomSelection{D} <: AbstractSelection
    "Distribution that is sampled."
    distribution::D
end

function (select::OrderedSelection)(prob, i, repeat)
    u0 = Dynamics.select_u0(prob.p, select.distribution, i)
    remake(prob, u0=u0)
end

function (select::RandomSelection)(prob, i, repeat)
    i = rand(1:length(select.distribution))
    u0 = Dynamics.select_u0(prob.p, select.distribution, i)
    remake(prob, u0=u0)
end
