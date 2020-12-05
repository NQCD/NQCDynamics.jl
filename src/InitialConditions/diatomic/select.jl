SUBROUTINE SELECT
# C         NSELT=2  CHOOSE INITIAL CONDITIONS FOR ONE OR TWO MOLECULES
# C         IN THE MAIN PROGRAM THE COORDINATES Q ARE SET EQUAL TO QZ
# C         AND THE MOMENTA P ARE SET EQUAL TO ZERO (FOR NSELT=2).
# think nthta is chosen mj but we dont choose it

#-->IN THIS CASE, ONLY REACTANT A IS MOVING, SURFACE IS FIXED AS REACTANT B 
system.shift_up()
#  ENERGY REFERENCE FOR SEPARATED REACTANTS
H,T,V = system.calc_energy()

V_ref=V

# C          DIATOM A IS TREATED SEMICLASSICALLY
#CALCULATE TURNING POINTS (done in the semiclassical quantization procedure.)

ENJA=(νᵢ+0.5)*DUM
RMIN,RMAX, ENJA = INITEBK(νᵢ,Jᵢ,RMIN,RMAX,V_ref,μ,ENJA,PTESTA,ALA)  #Gives RMIN and RMAX
SDUM=ENJA

# SELECT INITIAL RELATIVE COORDINATE AND MOMENTUM.

   #1. Generate random number
converged = false
do while (converged=false)
   RAND = gen(rand)
   #2. Calculate bond length, R
   R=RMIN+(RMAX-RMIN)*RAND
   #3. Calculate V(r)
   H,T,V = system.calc_energy(diatomic,bondlength=R)
   SUMM=ENJA-DUM/R**2-VDUM
   if (SUMM<=0.0)
      SUMM=0.0
      PR=0.0
   else
      #4. Calc p(r) , test if bigger than second random number, if not restart
      RAND = gen(rand)
      PR = sqrt(2*μ) * sqrt(ENJ - J(J+1)/(2μR**2) - V)
      SDUM=PTESTA/PR
      if (SDUM>RAND)
         converged = true
      end
   end
end
if (RAND<0.5) 
   PR=-PR
end
# CHOOSE INITIAL CARTESIAN COORDINATES AND MOMENTA, AND
# ANGULAR MOMENTUM.  DIATOM LIES ALONG THE X-AXIS.
# THEN RANDOMLY ROTATE THE CARTESIAN COORDINATES AND
# MOMENTA IN THE CENTER OF MASS FRAME.
CALL HOMOQP(R,PR,ALA,AM,μ,AI)

#CALCULATE ANGULAR MOMENTUM, MOMENT OF INERTIA TENSOR,
# ANGULAR VELOCITY, AND ROTATIONAL ENERGY
CALL ROTN(AM,EROTA)

function INITEBK(νᵢ,Jᵢ,RMIN,RMAX,V_ref,μ,ENJ,PTEST,AL)
   #    INITIALIZE PARAMETERS FOR AN OSCILLATOR WITH GIVEN
   #    QUANTUM NUMBERS N AND J BY SEMICLASSICAL EBK QUANTIZATION

   HNU=ENJ/(νᵢ+0.5)
   AM=√((Jᵢ*(Jᵢ+1)))
   AL=AM*ħ
   μ=calc_reduced_mass()
   DUM=1.0
   #  SOLVE FOR ENJ BY FIXED POINT APPROACH
   # i.e keep adjusting till you get vi requested
   i=0
   do while (abs(DUM)>.1.0D-6)
      RMIN,RMAX,AN,AJ = calc_NJ(ENJ,AM,RMIN,RMAX,V_ref,AN,AJ)
      DUM=νᵢ-AN
      ENJ=ENJ+DUM*HNU
      i=i+1
      if (i>1000)
        break
  end

   RMIN=RMIN+0.001
   RMAX=RMAX-0.001
   PTEST=√(0.0001*2.0*μ*ENJ)

   return
end


function HOMOQP(R,PR,AL,AM,μ,A)

   # SET CARTESIAN COORDINATES AND MOMENTA
   #Q positions, P momenta, R bondlength ?, μ reduced mass

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

   Q(J1)=R
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
   A(2)=μ*R**2
   A(3)=A(2)
   DUM=2π*RAND
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


function calc_NJ(ENJ,AM,RMIN,RMAX,V_ref,AN,AJ)
   using FastGaussQuadrature
   #    CALCULATE VIBRATIONAL AND ROTATIONAL QUANTUM NUMBERS FOR
   #    A PRODUCT DIATOM

   #     FIND THE ROTATIONAL QUANTUM NUMBER AJ FROM THE ROTATIONAL
   #     ANGULAR MOMENTUM AM IN UNITS OF H-BAR.

   AN=0
   AJ=0.5*(-1+√(1+4*AM²))
   #    FIND THE VIBRATIONAL QUANTUM NUMBER AN, USING AN INVERSION
   #    OF THE SEMICLASSICAL RYDBERG-KLEIN-REES (RKR) APPROACH.
   #
   #    INITIALIZE SOME VARIABLES FOR THE DIATOM.

   μ=calc_reduced_mass()
   T3=Q(2,3)-Q(1,3)
   T2=Q(2,2)-Q(1,2)
   T1=Q(2,1)-Q(1,1)
   RO=√(T1*T1+T2*T2+T3*T3)

   do i=1:3
      QO(i)=Q(3*L1-3+i)
      QO(i+3)=Q(3*L2-3+i)
      Q(3*L1-3+i)=(W(L1)*QO(i)+W(L2)*QO(i+3))/(W(L1)+W(L2))
      Q(3*L2-3+i)=(W(L1)*QO(i)+W(L2)*QO(i+3))/(W(L1)+W(L2))
   end
   Q(3*L1)=Q(3*L1)-0.5*RO
   Q(3*L2)=Q(3*L2)+0.5*RO
   H,T,V = diatomic.calc_energy(bondlength=R0)
   VEFF=V-V_ref+0.5*AM²*ħ²/(μ*RO²)

   # DETERMINE BOUNDARIES OF THE SEMICLASSICAL INTEGRAL
   do while (VEFF<ENJ and RZ<50.0)
      Q(3*L2)=Q(3*L2)+0.001
      H,T,V = system.calc_energy()
      RZ=Q(3*L2)-Q(3*L1)
      VEFF=V-V_ref+0.5*AM²*ħ²/(μ*RZ²)


   RMAX=RZ
   Q(3*L2)=Q(3*L1)+RO

   do while (VEFF<ENJ)
      Q(3*L2)=Q(3*L2)-0.001
      H,T,V = system.calc_energy()
      RZ=Q(3*L2)-Q(3*L1)
      VEFF=V-V_ref+0.5*AM²*ħ²/(μ*RZ²)

   RMIN=RZ

   # EVALUATE THE SEMICLASSICAL INTEGRAL BY GAUSSIAN QUADRATURE

   #Gauss-Legendre quadrature Parameters
   n_nodes = 50
   nodes, weights = gausslegendre(n_nodes)
   ASUM=0.0
   ΔR=(RMAX-RMIN)/50
   do i=1,n_nodes
      RZ=ΔR*(i-1)
      Q(3*L2)=Q(3*L1)+RZ
      H,T,V = system.calc_energy(Q)
      VEFF=(V-V_ref)+0.5*AM²*ħ²/(μ*RZ²)
      if (ENJ>VEFF)
         ASUM=ASUM+weights[i]*√(ENJ-VEFF)
   end
   do I=1,3
      Q(3*L1-3+I)=QO(I)
      Q(3*L2-3+I)=QO(I+3)
   end
   AN=√(8.0*μ)*ASUM/2π/ħ
   AN=AN-0.5

   end

