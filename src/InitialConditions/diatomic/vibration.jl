#Set momenta from selected initial vibrational state


function HOMOQP(R,PR,AL,AM,RMASS,A)

   # SET CARTESIAN COORDINATES AND MOMENTA
   #Q positions, P momenta, R bondlength ?, RMASS reduced mass?

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
   VREL=PR/RMASS
   WT=W(L(1))+W(L(2))
   VELA=VREL*W(L(2))/WT
   VELB=VELA-VREL
   J1=3*L(1)-2
   P(J1)=W(L(1))*VELA
   J1=3*L(2)-2
   P(J1)=W(L(2))*VELB
   # SET INERTIA ARRAYS.  CHOOSE Y AND Z ANGULAR MOMENTUM COMPONENTS
   A(1)=1.0D+20
   A(2)=RMASS*R**2
   A(3)=A(2)
   RAND=RAND0(ISEED)
   DUM=TWOPI*RAND
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
   WX=0.0D0
   WY=-AM(2)/A(2)
   WZ=-AM(3)/A(3)
   call ANGVEL(2)


   if (NTHTA.GE.0)
      if (NTHTA>JA)
            STOP 'NTHTA (Mj) CAN NOT EXCEED JA'
      else
         #  ROTATE A DIATOM (N=2) TO MODEL A SPECIFIC (J,M) STATE BY A VECTOR MODEL
         call ROTATEJM(NTHTA,JA,2)
      end
   else
      # RANDOMLY ROTATE THE DIATOM ABOUT ITS CENTER OF MASS BY
      # EULER'S ANGLES.  CENTER OF MASS COORDINATES QQ AND MOMENTA
      # PP ARE PASSED FROM SUBROUTINES CENMAS AND ANGVEL THROUGH
      # COMMON BLOCK WASTE.
      call ROTATE(2)
   end

   # ALCULATE ANGULAR MOMENTUM AND COMPONENTS.
   call ROTN(AM,EROT,2)
end

function INITEBK(N,J,RMIN,RMAX,DH,μ,ENJ,PTEST,AL)
    #    INITIALIZE PARAMETERS FOR AN OSCILLATOR WITH GIVEN
    #    QUANTUM NUMBERS N AND J BY SEMICLASSICAL EBK QUANTIZATION

    HNU=ENJ/(N+0.5)
    AM=√((J*(J+1)))
    AL=AM*ħ
    μ=calc_reduced_mass()
    DUM=1.0
    #  SOLVE FOR ENJ BY FIXED POINT APPROACH
    i=0
    do while (abs(DUM)>.1.0D-6)
       CALL calc_NJ(ENJ,AM,RMIN,RMAX,DH,AN,AJ)
       DUM=N-AN
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

function calc_NJ(ENJ,AM,RMIN,RMAX,DH,AN,AJ)

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
   L1=MIN(L(1),L(2))
   L2=MAX(L(1),L(2))
   μ=calc_reduced_mass()
   T3=Q(3*L2)-Q(3*L1)
   T2=Q(3*L2-1)-Q(3*L1-1)
   T1=Q(3*L2-2)-Q(3*L1-2)
   RO=√(T1*T1+T2*T2+T3*T3)

   do i=1:3
      QO(I)=Q(3*L1-3+I)
      QO(I+3)=Q(3*L2-3+I)
      Q(3*L1-3+I)=(W(L1)*QO(I)+W(L2)*QO(I+3))/(W(L1)+W(L2))
      Q(3*L2-3+I)=(W(L1)*QO(I)+W(L2)*QO(I+3))/(W(L1)+W(L2))
   end
   Q(3*L1)=Q(3*L1)-0.5*RO
   Q(3*L2)=Q(3*L2)+0.5*RO
   CALL energy
   VEFF=V-DH+0.5*AM²*ħ²/(μ*RO²)

   # DETERMINE BOUNDARIES OF THE SEMICLASSICAL INTEGRAL
   do while (VEFF<ENJ and RZ<50.0)
      Q(3*L2)=Q(3*L2)+0.001
      CALL energy
      RZ=Q(3*L2)-Q(3*L1)
      VEFF=V-DH+0.5*AM²*ħ²/(μ*RZ²)


   RMAX=RZ
   Q(3*L2)=Q(3*L1)+RO

   do while (VEFF<ENJ)
      Q(3*L2)=Q(3*L2)-0.001
      CALL energy
      RZ=Q(3*L2)-Q(3*L1)
      VEFF=V-DH+0.5*AM²*ħ²/(μ*RZ²)

   RMIN=RZ

   # EVALUATE THE SEMICLASSICAL INTEGRAL BY GAUSSIAN QUADRATURE
   NGL=50
   #GLPAR calculate parameters for Gauss-Legendre quadrature
   CALL GLPAR(RMIN,RMAX,XGL,WGL,NGL)
   ASUM=0.0
   do J=1,NGL
      RZ=XGL(J)
      Q(3*L2)=Q(3*L1)+RZ
      CALL ENERGY
      VEFF=V-DH+0.5*AM²*ħ²/(μ*RZ²)
      if (ENJ>VEFF)
         ASUM=ASUM+WGL(J)*√(ENJ-VEFF)
   end
   do I=1,3
      Q(3*L1-3+I)=QO(I)
      Q(3*L2-3+I)=QO(I+3)
   end
   AN=√(8.0*μ)*ASUM/2π/ħ
   AN=AN-0.5

end



