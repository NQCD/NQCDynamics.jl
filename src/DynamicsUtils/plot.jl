using RecipesBase: RecipesBase
using Dictionaries: Dictionary

RecipesBase.@recipe function f(dict::Dictionary, quantity::Symbol)
    xguide --> "t"
    yguide --> String(quantity)

    for i in eachindex(dict[quantity][1])
        RecipesBase.@series begin
            legend --> :false
            dict[:Time], [value[i] for value in table[quantity]]
        end
    end
end
