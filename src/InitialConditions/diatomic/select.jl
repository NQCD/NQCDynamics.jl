
#MAIN ROUTINE
#USES EVERYTHING
function select(νᵢJᵢ,Eᵢ,)
   # think nthta is chosen mj but we dont choose it
   #f is vibrational frequency
   #  ENERGY REFERENCE FOR SEPARATED REACTANTS
   H,T,V_ref = diatomic_calc_energy!(bondlength=5000,height=9)

   # DIATOM A IS TREATED SEMICLASSICALLY
   #CALCULATE TURNING POINTS (done in the semiclassical quantization procedure.)
   #f is vibrational frequency
   f = calc_force_constant()
   
   Evib=(νᵢ+0.5)*f*ħ   #desired vibrational energy

   limits, Evib, AL, AM, ptest = initialize_ebk!(νᵢ,Jᵢ,V_ref,Evib)  #Gives rmin and rmax, +parameters for ebk

   μ = diatomic_reduced_mass(R,atoms)


   # SELECT INITIAL RELATIVE COORDINATE AND MOMENTUM.
   converged = false
   while converged=false
      #1. Generate random number
      random = rand(1)
      #2. Calculate bond length, r
      r=rmin+(rmax-rmin)*random
      #3. Calculate V(r)
      H,T,V = diatomic.calc_energy(bondlength=r)
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
   calc_ROTN!(AM,EROTA)

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
   return
end

