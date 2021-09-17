using RecipesBase: @recipe, @series
using TypedTables: Table

@recipe function f(table::Table, quantity::Symbol)
    xguide --> "t"
    yguide --> String(quantity)

    for i in eachindex(eval(:($table.$quantity[1])))
        @series begin
            legend --> :false
            table.t, eval(:([value[$i] for value in $table.$quantity]))
        end
    end
end
