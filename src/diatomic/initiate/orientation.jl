#Randomly rotate diatomic around COM, rotate momenta accordingly
#In future allow specific range of angles to be selected e.g 0-60 degrees

function orient(system,R::positions, P :: momenta)

    θ = cos^-1(2*rand1-1)
    ϕ = 2*π*rand2
    χ = 2*π*rand3

    COM,μ = calc_COM()
    COM_P = calc_COM_momenta()
    
    RXX=cos(θ)*cos(ϕ)*cos(χ)-sin(ϕ)*sin(χ)
    RXY=-cos(θ)*cos(ϕ)*sin(χ)-sin(ϕ)*cos(χ)
    RXZ=sin(θ)*cos(ϕ)
    RYX=cos(θ)*sin(ϕ)*cos(χ)+cos(ϕ)*sin(χ)
    RYY=-cos(θ)*sin(ϕ)*sin(χ)+cos(ϕ)*cos(χ)
    RYZ=sin(θ)*sin(ϕ)
    RZX=-sin(θ)*cos(χ)
    RZY=sin(θ)*sin(χ)
    RZZ=cos(θ)
    

    for i in 1:N
        J=3*L[i]
        R[J-2]=COM[J-2]*RXX+COM[J-1]*RXY+COM[J]*RXZ
        R[J-1]=COM[J-2]*RYX+COM[J-1]*RYY+COM[J]*RYZ
        R[J]=COM[J-2]*RZX+COM[J-1]*RZY+COM[J]*RZZ
        COM[J-2]=R[J-2]
        COM[J-1]=R[J-1]
        COM[J]=R[J] 
        P[J-2]=COM_P[J-2]*RXX+COM_P[J-1]*RXY+COM_P[J]*RXZ
        P[J-1]=COM_P[J-2]*RYX+COM_P[J-1]*RYY+COM_P[J]*RYZ
        P[J]=COM_P[J-2]*RZX+COM_P[J-1]*RZY+COM_P[J]*RZZ
        COM_P[J-2]=P[J-2]
        COM_P[J-1]=P[J-1]
        COM_P[J]=P[J] 
    end
end