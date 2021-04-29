using DiffEqBase: CallbackSet

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
    v, r = InitialConditions.pick(select.distribution, i)
    u0 = select_u0(prob.p, v, r, select.distribution.state)
    remake(prob, u0=u0)
end

function (select::RandomSelection)(prob, i, repeat)
    v, r = rand(select.distribution)
    u0 = select_u0(prob.p, v, r, select.distribution.state)
    remake(prob, u0=u0)
end

struct SelectWithCallbacks{S<:AbstractSelection,C1,C2,V} <: AbstractSelection
    selection::S
    standard_callbacks::C1
    changing_callbacks::C2
    values::V
    function SelectWithCallbacks(selection, standard_callbacks, output, trajectories)

        callbacks = []
        values = []
        for i=1:trajectories
            cb, vals = Dynamics.create_saving_callback(output)
            push!(callbacks, cb)
            push!(values, vals)
        end

        new{typeof(selection), typeof(standard_callbacks), typeof(callbacks), typeof(values)}(
            selection, standard_callbacks, callbacks, values)
    end
end

function (select::SelectWithCallbacks)(prob, i, repeat)
    prob = select.selection(prob, i, repeat)
    remake(prob, callback=CallbackSet(select.standard_callbacks, select.changing_callbacks[i]))
end
