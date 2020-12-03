#Set momenta from selected initial vibrational state

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
   H,T,V = system.calc_energy()
   VEFF=V-DH+0.5*AM²*ħ²/(μ*RO²)

   # DETERMINE BOUNDARIES OF THE SEMICLASSICAL INTEGRAL
   do while (VEFF<ENJ and RZ<50.0)
      Q(3*L2)=Q(3*L2)+0.001
      H,T,V = system.calc_energy()
      RZ=Q(3*L2)-Q(3*L1)
      VEFF=V-DH+0.5*AM²*ħ²/(μ*RZ²)


   RMAX=RZ
   Q(3*L2)=Q(3*L1)+RO

   do while (VEFF<ENJ)
      Q(3*L2)=Q(3*L2)-0.001
      H,T,V = system.calc_energy()
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
      H,T,V = system.calc_energy()
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



