function select(νᵢJᵢ,Eᵢ,)
   # think nthta is chosen mj but we dont choose it
   #f is vibrational frequency
   #  ENERGY REFERENCE FOR SEPARATED REACTANTS
   H,T,V_ref = diatomic_calc_energy!(bondlength=5000,height=9)

   # DIATOM A IS TREATED SEMICLASSICALLY
   #CALCULATE TURNING POINTS (done in the semiclassical quantization procedure.)
   enj=(νᵢ+0.5)*f*ħ   #desired vibrational energy
   μ = calc_reduced_mass()
   limits, enj, AL, AM, ptest = INITEBK!(νᵢ,Jᵢ,V_ref,enj)  #Gives rmin and rmax

   # SELECT INITIAL RELATIVE COORDINATE AND MOMENTUM.

      #1. Generate random number
   converged = false
   while converged=false
      random = rand(1)
      #2. Calculate bond length, r
      r=rmin+(rmax-rmin)*random
      #3. Calculate V(r)
      H,T,V = system.calc_energy(diatomic,bondlength=r)
      SUMM=enj-DUM/r^2-VDUM
      if (SUMM<=0.0)
         SUMM=0.0
         PR=0.0
      else
         #4. Calc p(r) , test if bigger than second random number, if not restart
         random = rand(1)
         PR = sqrt(2*μ) * sqrt(enj - J(J+1)/(2μR^2) - V)
         SDUM=ptest/PR
         if SDUM>random
            converged = true
         end
      end
   end
   if (random<0.5) 
      PR=-PR
   end
   # CHOOSE INITIAL CARTESIAN COORDINATES AND MOMENTA, AND
   # ANGULAR MOMENTUM.  DIATOM LIES ALONG THE X-AXIS.
   # THEN RANDOMLY ROTATE THE CARTESIAN COORDINATES AND
   # MOMENTA IN THE CENTER OF MASS FRAME.
   HOMOQP!(r,PR,AL,AM,μ,AI)

   #CALCULATE ANGULAR MOMENTUM, MOMENT OF INERTIA TENSOR,
   # ANGULAR VELOCITY, AND ROTATIONAL ENERGY
   ROTN!(AM,EROTA)

end

function INITEBK!(νᵢ,Jᵢ,V_ref,enj)
   #    INITIALIZE PARAMETERS FOR AN OSCILLATOR WITH GIVEN
   #    QUANTUM NUMBERS n AND J BY SEMICLASSICAL EBK QUANTIZATION

   μ = diatomic_reduced_mass()
   HNU=enj/(νᵢ+0.5)
   AM=√((Jᵢ*(Jᵢ+1)))
   AL=AM*ħ
   
   #  SOLVE FOR enj BY FIXED POINT APPROACH
   # i.e keep adjusting till you get vi requested
   DUM=1.0
   i=0
   while (abs(DUM)>1.0e-6)
       n,J,limits = calc_NJ(enj,AM,V_ref)
       DUM=νᵢ-n
       enj=enj+DUM*HNU
       i=i+1
       if (i>1000)
           break
       end
   end


   limits[1]=limits[1]+0.001
   limits[2]=limits[2]-0.001
   ptest=√(0.0001*2.0*μ*enj)

   return limits, enj, AL, AM, ptest
end


function HOMOQP!(r,PR,AL,AM,μ,A)

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
   call CENMAS(WT,QCM,VCM,2)
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
      call ROTATE(2)

   # CALCULATE ANGULAR MOMENTUM AND COMPONENTS.
   call ROTN(AM,EROT,2)
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

#Sets translational momenta for COM
#For now this is fixed input parameter
#currently a mess not sure what im doing

function translate(system,positions,v,ei,θᵢ)
   #ei is user defined incidence translation energy
   #theta_i is user defined incidence angle
   masses=system.get_masses()

   com,μ = calc_COM(system)  

   com_v=[0,0,-sqrt(2*ei/μ)] #straight down

   com_v=[0,com_v[2]*cos(θᵢ)-com_v[3]*sin(θᵢ),com_v[2]*sin(θᵢ)-com_v[3]*cos(θᵢ)] #rotate by θᵢ
   
   v[1,:] += com_v
   v[2,:] += com_v

   return v
end

 #Displace diatomic to specified height
 
function height(system,R :: positions,H :: height)
   #todo: allow any surface plane by allowing user defining orthogonal vector to surface plane
   COM,μ = calc_COM(system)

   COM[2] = height_original

   displacement = height - height_original

   positions(2,:) += displacement

end

