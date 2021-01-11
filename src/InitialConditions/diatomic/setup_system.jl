#USES calc_com, ANGVEL, orient, calc_ROTN
#NEEDS WORK
function set_positions_momenta!(r,PR,AL,AM,μ,A)
    #Prior name: HOMOQP
    # SET CARTESIAN COORDINATES AND MOMENTA
    #Q positions, P momenta, r bondlength ?, μ reduced mass
 
    for I=1:2
       J3=3*L(I)
       J2=J3-1
       J1=J2-1
       Q(J1)=0.0
       Q(J2)=0.0
       Q(J3)=0.0
       P(J2)=0.0
       P(J3)=0.0
    end
 
    Q(J1)=r
    VREL=PR/μ
    WT=W(L(1))+W(L(2))
    VELA=VREL*W(L(2))/WT
    VELB=VELA-VREL
    J1=3*L(1)-2
    P(J1)=W(L(1))*VELA
    J1=3*L(2)-2
    P(J1)=W(L(2))*VELB
 
 
    # SET INERTIA ARRAYS.  CHOOSE Y AND Z ANGULAR MOMENTUM COMPONENTS
    A(1)=1.0D+20
    A(2)=μ*r**2
    A(3)=A(2)
    DUM=2π*random
    AM(1)=0.0D0
    AM(2)=AL*SIN(DUM)
    AM(3)=AL*COS(DUM)
    if (NTHTA>0) THEN
       AM(3)=AL
       AM(2)=0.0D0
    end
    #  CALCULATE CENTER OF MASS COORDIANTES QQ AND MOMENTA PP
    COM_R, COM_P = calculate_COM(WT,QCM,VCM,2)
    #  MOVE PP ARRAY TO P ARRAY AND QQ ARRAY TO Q ARRAY
    for I=1:2
       J=3*L(I)+1
       for K=1:3
          Q(J-K)=QQ(J-K)
          P(J-K)=PP(J-K)
       end
    end
    # ADD THE ANGULAR MOMENTUM VECTOR.  CALCULATE THE REQUIRED
    # ANGULAR VELOCITY AND ADD IT TO THE DIATOM, WHICH LIES
    # ALONG THE X-AXIS.
    WX=0.0
    WY=-AM(2)/A(2)
    WZ=-AM(3)/A(3)
    call ANGVEL(2)
 
    # RANDOMLY ROTATE THE DIATOM ABOUT ITS CENTER OF MASS BY
    # EULER'S ANGLES.  CENTER OF MASS COORDINATES QQ AND MOMENTA
    # PP ARE PASSED FROM SUBROUTINES CENMAS AND ANGVEL THROUGH
    # COMMON BLOCK WASTE.
    call orient(2)
 
    # CALCULATE ANGULAR MOMENTUM AND COMPONENTS.
    call calc_ROTN(AM,EROT,2)
 end

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
   
 
   for i in 1:2
       J=3*L[i] #what is L?
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