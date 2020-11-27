#Set momenta from selected initial vibrational state

function set_vibration(system,positions)

    COM,μ = calc_COM(system,positions)

    ρ1 = +(2*μ )^(1/2) * (E - J(J+1)/2*μr^2 - V )^(1/2)
    ρ2 = -1*p1

    return positions,momenta