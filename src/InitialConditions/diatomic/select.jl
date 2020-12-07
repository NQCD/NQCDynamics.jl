SUBROUTINE SELECT
# C         NSELT=2  CHOOSE INITIAL CONDITIONS FOR ONE OR TWO MOLECULES
# C         IN THE MAIN PROGRAM THE COORDINATES Q ARE SET EQUAL TO QZ
# C         AND THE MOMENTA P ARE SET EQUAL TO ZERO (FOR NSELT=2).
# think nthta is chosen mj but we dont choose it
#f is vibrational frequency
#  ENERGY REFERENCE FOR SEPARATED REACTANTS
H,T,V_ref = diatomic_calc_energy!(bondlength=5000,height=9)

# C          DIATOM A IS TREATED SEMICLASSICALLY
#CALCULATE TURNING POINTS (done in the semiclassical quantization procedure.)

enj=(νᵢ+0.5)*f*ħ   #desired vibrational energy
μ = calc_reduced_mass()
rmin,rmax, enj = INITEBK!(νᵢ,Jᵢ,rmin,rmax,V_ref,μ,enj,ptest,ALA)  #Gives rmin and rmax

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
HOMOQP!(r,PR,ALA,AM,μ,AI)

#CALCULATE ANGULAR MOMENTUM, MOMENT OF INERTIA TENSOR,
# ANGULAR VELOCITY, AND ROTATIONAL ENERGY
ROTN!(AM,EROTA)

function INITEBK!(νᵢ,Jᵢ,rmin,rmax,V_ref,μ,enj,ptest,AL)
   #    INITIALIZE PARAMETERS FOR AN OSCILLATOR WITH GIVEN
   #    QUANTUM NUMBERS n AND J BY SEMICLASSICAL EBK QUANTIZATION

   HNU=enj/(νᵢ+0.5)
   AM=√((Jᵢ*(Jᵢ+1)))
   AL=AM*ħ
   
   #  SOLVE FOR enj BY FIXED POINT APPROACH
   # i.e keep adjusting till you get vi requested
   DUM=1.0
   i=0
   while (abs(DUM)>1.0e-6)
       rmin,rmax,n,J = calc_NJ(enj,AM,rmin,rmax,V_ref,n,J)
       DUM=νᵢ-n
       enj=enj+DUM*HNU
       i=i+1
       if (i>1000)
           break
       end
   end

   rmin=rmin+0.001
   rmax=rmax-0.001
   ptest=√(0.0001*2.0*μ*enj)

   return nothing
end


function HOMOQP!(r,PR,AL,AM,μ,A)

   # SET CARTESIAN COORDINATES AND MOMENTA
   #Q positions, P momenta, r bondlength ?, μ reduced mass

   do I=1:2
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
   do I=1,2
      J=3*L(I)+1
      do K=1,3
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


function calc_NJ!(enj,μ,AM,rmin,rmax,V_ref,n,J)
   using FastGaussQuadrature
   #    CALCULATE VIBRATIONAL AND ROTATIONAL QUANTUM NUMBERS FOR
   #    A PRODUCT DIATOM

   #     FIND THE ROTATIONAL QUANTUM NUMBER J FROM THE ROTATIONAL
   #     ANGULAR MOMENTUM AM IN UNITS OF ħ.

   n=0
   J=0.5*(-1+√(1+4*AM²))
   #    FIND THE VIBRATIONAL QUANTUM NUMBER n, USING AN INVERSION
   #    OF THE SEMICLASSICAL RYDBERG-KLEIN-REES (RKR) APPROACH.
   #
   #    INITIALIZE SOME VARIABLES FOR THE DIATOM.

   T3=Q(2,3)-Q(1,3)
   T2=Q(2,2)-Q(1,2)
   T1=Q(2,1)-Q(1,1)
   RO=√(T1*T1+T2*T2+T3*T3)

   H,T,V = diatomic.calc_energy(bondlength=R0)
   VEFF=V-V_ref+0.5*AM²*ħ²/(μ*RO²)

   # DETERMINE BOUNDARIES OF THE SEMICLASSICAL INTEGRAL
   while VEFF<enj & RZ<50.0
       RZ=RZ+0.001
       H,T,V = diatomic.calc_energy(bondlength=RZ)
       VEFF=V-V_ref+0.5*AM²*ħ²/(μ*RZ²)
   end
   rmax=RZ

   RZ=R0
   while VEFF<enj
      RZ=RZ-0.001
      H,T,V = diatomc.calc_energy(bondlength=RZ)
      VEFF=V-V_ref+0.5*AM²*ħ²/(μ*RZ²)
   end
   rmin=RZ

   # EVALUATE THE SEMICLASSICAL INTEGRAL BY GAUSSIAN QUADRATURE

   #Gauss-Legendre quadrature Parameters
   n_nodes = 50
   nodes, weights = gausslegendre(n_nodes)
   ASUM=0.0
   ΔR=(rmax-rmin)/50
   for i=1:n_nodes
      RZ=ΔR*(i-1)
      R0 = rmin+RZ
      H,T,V = diatomic.calc_energy(bondlength=R0)
      VEFF=(V-V_ref)+0.5*AM²*ħ²/(μ*RZ²)
      if (enj>VEFF)
         ASUM=ASUM+weights[i]*√(enj-VEFF)
      end
   end

   n=√(8.0*μ)*ASUM/2π/ħ
   n=n-0.5

end