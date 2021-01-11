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
        n,J,limits = calc_NJ(Evib,AM,V_ref)
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