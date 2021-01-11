
#MAIN ROUTINE
#USES EVERYTHING
function select(νᵢJᵢ,Eᵢ,)
   # think nthta is chosen mj but we dont choose it
   #f is vibrational frequency
   #  ENERGY REFERENCE FOR SEPARATED REACTANTS
   H,T,V_ref = calculate_diatomic_energy!(bondlength=5000,height=9)

   # DIATOM A IS TREATED SEMICLASSICALLY
   #CALCULATE TURNING POINTS (done in the semiclassical quantization procedure.)
   #f is vibrational frequency
   f = calculate_force_constant()
   
   Evib=(νᵢ+0.5)*f*ħ   #desired vibrational energy

   limits, Evib, AL, AM, ptest = initialize_ebk!(νᵢ,Jᵢ,V_ref,Evib)  #Gives rmin and rmax, +parameters for ebk

   μ = calculate_reduced_mass(R,atoms)


   # SELECT INITIAL RELATIVE COORDINATE AND MOMENTUM.
   converged = false
   while converged=false
      #1. Generate random number
      random = rand(1)
      #2. Calculate bond length, r
      r=rmin+(rmax-rmin)*random
      #3. Calculate V(r)
      H,T,V = calculate_diatomic_energy(bondlength=r)
      SUMM=Evib-DUM/r^2-VDUM
      if (SUMM<=0.0)
         SUMM=0.0
         PR=0.0
      else
         #4. Calc p(r) , test if bigger than second random number, if not restart
         random = rand(1)
         PR = sqrt(2*μ) * sqrt(Evib - J(J+1)/(2μR^2) - V)
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
   set_positions_momenta!(r,PR,AL,AM,μ,AI)

   #CALCULATE ANGULAR MOMENTUM, MOMENT OF INERTIA TENSOR,
   # ANGULAR VELOCITY, AND ROTATIONAL ENERGY
   calculate_ROTN!(AM,EROTA)

end

#USES calc_NJ
#BASICALLY READY FOR TESTING
function initialize_ebk!(νᵢ,Jᵢ,V_ref,Evib)
   #Prior name: INITEBK
   #    INITIALIZE PARAMETERS FOR AN OSCILLATOR WITH GIVEN
   #    QUANTUM NUMBERS n AND J BY SEMICLASSICAL EBK QUANTIZATION

   μ = calc_μ(atoms)
   HNU=Evib/(νᵢ+0.5)
   AM=sqrt((Jᵢ*(Jᵢ+1)))
   AL=AM*ħ
   
   #  SOLVE FOR Evib BY FIXED POINT APPROACH
   # i.e keep adjusting till you get vi requested
   diff=1.0
   tolerance = 1.0e-6
   i=0 #iter counter
   while (abs(diff)>tolerance) 
       n,J,limits = calculate_NJ(Evib,AM,V_ref)
       diff=νᵢ-n
       Evib=Evib+diff*HNU
       i=i+1
       if (i>1000) #max iter
           break
       end
   end


   limits[1]=limits[1]+0.001
   limits[2]=limits[2]-0.001
   ptest=sqrt(0.0001*2.0*μ*Evib)

   limits, Evib, AL, AM, ptest
end

function ANGVEL(N)
   # C         SUBTRACT OFF THE ANGULAR VELOCITY
   # C
   do I=1,N
      J=L(I)
      J3=3*J
      J2=J3-1
      J1=J2-1
      P(J1)=P(J1)-(QQ(J3)*WY-QQ(J2)*WZ)*W(J)
      P(J2)=P(J2)-(QQ(J1)*WZ-QQ(J3)*WX)*W(J)
      P(J3)=P(J3)-(QQ(J2)*WX-QQ(J1)*WY)*W(J)
      PP(J1)=P(J1)
      PP(J2)=P(J2)
      PP(J3)=P(J3)
   end
end

