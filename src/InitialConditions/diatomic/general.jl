function calc_COM(R :: positions)
    #calculate centre of mass from diatomic coordinates


    return COM
end

function calc_COM_momenta(P :: momenta)

    return COM_P
end

function calc_diatomic_energy(bondlength=R)
    # wrapper function to calc diatomic potential and kinetic energy
    # with bondlength R and height from surface Z
    height = 9 #large enough to minimize interaction, small enough to be in training
                #space

    return H, T, V
end

function calc_reduced_mass()
    #calc reduced mass (μ)


    return μ
end