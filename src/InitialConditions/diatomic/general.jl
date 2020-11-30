function calc_COM(R :: positions)
    #calculate centre of mass from diatomic coordinates


    return COM
end

function calc_COM_momenta(P :: momenta)

    return COM_P
end

function energy(H,T,V)
    # CALCULATE POTENTIAL, KINETIC AND TOTAL ENERGY OF THE
    # MOLECULAR SYSTEM

    T=0.0
    V=0.0
  
    CALL POT0(NATOMS,VV)
    #  ADD VZERO TO THE POTENTIAL ENERGY
    V=VV+VZERO
    #        CALCULATE KINETIC ENERGY
    J=1
    DO I=1,NATOMS
    T=T+(P(J)**2+P(J+1)**2+P(J+2)**2)/2.0/W(I)
    J=J+3
    ENDDO

    H=T+V
end

function calc_reduced_mass()
    #calc reduced mass (μ)


    return μ
end