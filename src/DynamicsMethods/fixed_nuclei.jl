using Interpolations
# For diabatic methods it seems that the electornics are updated via set_quantum_derivative and the update_electornics! function
# Should be able to write a motion! loop that only calls these functions
# Not too sure about mapping variables so will ignore that for now
# This appears to still be the same workflow for ring-polymer simulations - not that I see a reason to do frozne nuclei ring polymer

function DynamicsMethods.frozennuclei!(du, u, p::Tuple{Real,Simulation}, t)
    r = p[1]
    sim = p[2]
    Calcultors.update_electronics!(sim.calculator, r)
    DynamicsUtils.set_quantum_derivative!(du, u, sim) # 
end

function DynamicsMethods.frozennuclei!(du, u, p::Tuple{Interpolations.AbstractInterpolation
    ,Simulation}, t)
    r = p[1](t)
    sim = p[2]
    Calcultors.update_electronics!(sim.calculator, r)
    DynamicsUtils.set_quantum_derivative!(du, u, sim)
end

# That should be everything required honestly. Now just need to write a keyword which sets nuclear position as a parameter rather than a variable and 
# calls frozennuclei! rather than motion! when constructing the problem
# Thanks for the input :) - I make notes all over what I code like this sorry
