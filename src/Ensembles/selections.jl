using DiffEqBase: CallbackSet

using NonadiabaticMolecularDynamics: InitialConditions, DynamicsOutputs

abstract type AbstractSelection end

"""
Select the initial conditions from the distribution in order. 
"""
struct OrderedSelection{D,I} <: AbstractSelection
    "Distribution that is sampled."
    distribution::D
    indices::I
end

"""
Obtain initial conditions by randomly sampling the distribution.
"""
struct RandomSelection{D} <: AbstractSelection
    "Distribution that is sampled."
    distribution::D
end

function (select::OrderedSelection)(prob, i, repeat)
    j = select.indices[i]
    u = InitialConditions.pick(select.distribution, j)
    u0 = select_u0(prob.p, u.v, u.r, select.distribution.state, select.distribution.type)
    DynamicsMethods.create_problem(u0, prob.tspan, prob.p)
end

function (select::RandomSelection)(prob, i, repeat)
    u = rand(select.distribution)
    u0 = select_u0(prob.p, u.v, u.r, select.distribution.state, select.distribution.type)
    DynamicsMethods.create_problem(u0, prob.tspan, prob.p)
end

struct SelectWithCallbacks{S<:AbstractSelection,C1,C2,V} <: AbstractSelection
    selection::S
    standard_callbacks::C1
    changing_callbacks::C2
    values::V
    function SelectWithCallbacks(selection, standard_callbacks, output, trajectories; saveat=[])

        callbacks = []
        values = []
        for i=1:trajectories
            cb, vals = DynamicsOutputs.create_saving_callback(output; saveat=saveat)
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
