#Displace diatomic to specified height

function height(system,R :: positions,H :: height)
    #todo: allow any surface plane by allowing user defining orthogonal vector to surface plane
    COM,μ = calc_COM(system)

    COM[2] = height_original

    displacement = height - height_original

    positions(2,:) += displacement

end