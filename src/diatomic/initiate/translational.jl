#Sets translational momenta for COM
#For now this is fixed input parameter
#currently a mess not sure what im doing

function translate(system,positions,v,ei,θᵢ)
    #ei is user defined incidence translation energy
    #theta_i is user defined incidence angle
    masses=system.get_masses()

    com,rm = calc_COM(system)  

    com_v=[0,0,-sqrt(2*ei/rm)] #straight down

    com_v=[0,com_v[2]*cos(θᵢ)-com_v[3]*sin(θᵢ),com_v[2]*sin(θᵢ)-com_v[3]*cos(θᵢ)] #rotate by θᵢ
    
    v[1,:] += com_v
    v[2,:] += com_v

    return v
end