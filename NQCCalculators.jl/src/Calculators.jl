#Definitions of base NQCModels functions on an Abstract Cache.
NQCModels.nstates(cache::Abstract_Cache) = NQCModels.nstates(cache.model)
NQCModels.eachstate(cache::Abstract_Cache) = NQCModels.eachstate(cache.model)
NQCModels.nelectrons(cache::Abstract_Cache) = NQCModels.nelectrons(cache.model)
NQCModels.eachelectron(cache::Abstract_Cache) = NQCModels.eachelectron(cache.model)
NQCModels.mobileatoms(cache::Abstract_Cache) = NQCModels.mobileatoms(cache.model, size(cache.derivative, 2))
NQCModels.dofs(cache::Abstract_Cache) = NQCModels.dofs(cache.model)
NQCModels.fermilevel(cache::Abstract_Cache) = NQCModels.fermilevel(cache.model)
beads(cache) = Base.OneTo(length(cache.potential))
Base.eltype(::Abstract_Cache{T}) where {T} = T

"""
Each of the quantities specified here has functions:
`get_quantity(cache, r)`
`evaluate_quantity!(cache, r)`
`update_quantity!(cache, r)`

The user should mostly access the get_quantity!() function, or in rare circumstances the evaluate_quantity!() fucntion.
This will ensure quantities are correctly evaluated and cached accordingly.

The latter is called by the former and is where the details required to calculate the quantity are found.
"""
const quantities = [
    :potential,
    :derivative,
    :eigen,
    :adiabatic_derivative,
    :nonadiabatic_coupling,

    :traceless_potential,
    :V̄,
    :traceless_derivative,
    :D̄,
    :traceless_adiabatic_derivative,

    :centroid,
    :centroid_potential,
    :centroid_derivative,
    :centroid_eigen,
    :centroid_adiabatic_derivative,
    :centroid_nonadiabatic_coupling,

    :friction,
]

for quantity in quantities
    get_quantity = Symbol(:get_, quantity)
    field = Expr(:call, :getfield, :cache, QuoteNode(quantity))

    @eval function $(get_quantity)(cache, r)
        return $field
    end
end

include("Evaluate_Functions.jl")

include("Update_Functions.jl")