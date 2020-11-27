#Set bond length

function set_bond_length(system,positions)

    r1,r2 = find_turning_points(system,positions)
    r = r1 + (r2-r1)*rand1


    return positions
end


function find_turning_points(system,positions)
    #"semiclassical quantization procedure"

    return r1,r2
end