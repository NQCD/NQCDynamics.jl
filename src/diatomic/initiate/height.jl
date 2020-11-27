#Displace diatomic to specified height
#Not sure if this should be apart of it


function height(system,positions,height)

    COM,μ = calc_COM(system)

    COM[2] = height_original

    displacement = height - height_original

    positions(2,:) += displacement

end