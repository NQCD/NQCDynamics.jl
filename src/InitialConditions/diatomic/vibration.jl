#Set momenta from selected initial vibrational state

function set_vibration(N,J,RMIN,RMAX,DH,\mu,ENJ,PTEST,AL)
    #    INITIALIZE PARAMETERS FOR AN OSCILLATOR WITH GIVEN
    #    QUANTUM NUMBERS N AND J BY SEMICLASSICAL EBK QUANTIZATION

    HNU=ENJ/(N+0.5)
    AM=√((J*(J+1)))
    AL=AM*ħ
    μ=calc_reduced_mass()
    DUM=1.0
    #  SOLVE FOR ENJ BY FIXED POINT APPROACH
    ICOUNT=0
    DO WHILE (abs(DUM)>.1.0D-6)
       CALL calc_NJ(ENJ,AM,RMIN,RMAX,DH,AN,AJ)
       DUM=N-AN
       ENJ=ENJ+DUM*HNU
       ICOUNT=ICOUNT+1
       IF (ICOUNT>1000) then
         break
    END DO

    RMIN=RMIN+0.001
    RMAX=RMAX-0.001
    PTEST=√(0.0001*2.0*μ*ENJ)

    return
end

function calc_NJ(ENJ,AM,RMIN,RMAX,DH,AN,AJ)

    #         CALCULATE VIBRATIONAL AND ROTATIONAL QUANTUM NUMBERS FOR
    #         A PRODUCT DIATOM

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

   DO I=1,3
      QO(I)=Q(3*L1-3+I)
      QO(I+3)=Q(3*L2-3+I)
      Q(3*L1-3+I)=(W(L1)*QO(I)+W(L2)*QO(I+3))/(W(L1)+W(L2))
      Q(3*L2-3+I)=(W(L1)*QO(I)+W(L2)*QO(I+3))/(W(L1)+W(L2))
   ENDDO
   Q(3*L1)=Q(3*L1)-0.5*RO
   Q(3*L2)=Q(3*L2)+0.5*RO
   CALL energy
   VEFF=V-DH+0.5*AM²*ħ²/(μ*RO²)

   # DETERMINE BOUNDARIES OF THE SEMICLASSICAL INTEGRAL
   do while (VEFF<ENJ and RZ<50.0)
      Q(3*L2)=Q(3*L2)+0.001
      CALL energy
      RZ=Q(3*L2)-Q(3*L1)
      VEFF=V-DH+0.5*AM²*ħ*ħ/(C1*\mu*RZ*RZ)


   RMAX=RZ
   Q(3*L2)=Q(3*L1)+RO

   do while (VEFF<ENJ)
      Q(3*L2)=Q(3*L2)-0.001
      CALL energy
      RZ=Q(3*L2)-Q(3*L1)
      VEFF=V-DH+0.5*AM²*ħ*ħ/(μ*RZ*RZ)


   RMIN=RZ
   # EVALUATE THE SEMICLASSICAL INTEGRAL BY GAUSSIAN QUADRATURE
   NGL=50
   CALL GLPAR(RMIN,RMAX,XGL,WGL,NGL)
   ASUM=0.0
   do J=1,NGL
      RZ=XGL(J)
      Q(3*L2)=Q(3*L1)+RZ
      CALL ENERGY
      VEFF=V-DH+0.5*AM²*ħ²/(C1*μ*RZ*RZ)
      IF (ENJ>VEFF)
         ASUM=ASUM+WGL(J)*√(ENJ-VEFF)
   end
   do I=1,3
      Q(3*L1-3+I)=QO(I)
      Q(3*L2-3+I)=QO(I+3)
   end
   AN=√(8.0*μ)*ASUM/2π/ħ
   AN=AN-0.5

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

end