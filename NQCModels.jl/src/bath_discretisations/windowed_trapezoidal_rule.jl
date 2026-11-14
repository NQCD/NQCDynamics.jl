"""
    WindowedTrapezoidalRule{B,T}

Discretise wide band continuum using the trapezoidal rule with a region of higher state density
surrounding the Fermi Energy. The position of the energy window containing more densely packed states
is defined by energy values and the increase in density by a proportion. 
"""

struct WindowedTrapezoidalRule{T} <: WideBandBathDiscretisation
    bathstates::Vector{T} # due to the non-even distribution need to represent as vectors
    bathcoupling::Vector{T} 
end

"""
    WindowedTrapezoidalRule(M, bandmin, bandmax, windmin, windmax; densityratio=0.50, fermilevelstate=false)
               
"""
function WindowedTrapezoidalRule(M, bandmin, bandmax, windmin, windmax; densityratio=0.50, fermilevelstate=false)
    (M - (M*densityratio)) % 2 == 0 || throw(error("For the provided arguments, the number of states not in the windowed region is $(M - (M*densityratio)). This value must be even for a symmetric discretisation.")) # constraint enforced such that the a densityratio=0.5 can be utilized and return an integer number of states.
    abs(bandmin) > abs(windmin) || throw(error("Requested window minimum energy lies outside of energy range."))
    abs(bandmax) > abs(windmax) || throw(error("Requested window maximum energy lies outside of energy range."))
    
    M_window = floor(Int, M*densityratio)
    M_sparse = floor(Int, (M - (M*densityratio))/2)

    if fermilevelstate == false
        nothing
    elseif fermilevelstate == true
        M_window = M_window+1 # increase number of states in window region by 1 to place a state at the fermi level
        M_new = M_window + M_sparse * 2
        @info "Additional state added at Fermi level, number of states in discretisation =$M_new"
    else
        throw(error("Invalid value provided to `fermilevelstate`. Only Boolean values allowed."))
    end

    # densely distributed trapezoidal rule discretised bath states within energy window
    ΔE_window = windmax - windmin
    bstates_window = collect(range(windmin, windmax, length=M_window))

    # window spacing
    δE_window = bstates_window[2] - bstates_window[1]

    # sparse region energy spacing 
    δE_sparse = 2*(bandmax - windmax - 0.5*δE_window) / (2*M_sparse -1)
    
    # sparse to window energy spacing
    δE_join = 0.5*(δE_sparse + δE_window)

    # Sparsely distributed trapezoidal rule discretised bath states
    ΔE_sparse1 = windmin - δE_join  - bandmin # Energy range for first sparsely distributed state region
    bstates_sparse1 = collect(range(bandmin, windmin - δE_join, length=M_sparse)) 

    ΔE_sparse2 = bandmax - windmax - δE_join # Energy Range for second sparsely distributed state region
    bstates_sparse2 = collect(range(windmax + δE_join, bandmax, length=M_sparse))
    

    bathstates = [bstates_sparse1; bstates_window; bstates_sparse2] # concatenates arrays vertically (along axis = 1)    
    
    bcoupling_window = fill!(copy(bstates_window), sqrt(ΔE_window / (M_window-1)))
    bcoupling_sparse1 = fill!(copy(bstates_sparse1), sqrt(ΔE_sparse1 / (M_sparse-1)))
    bcoupling_sparse2 = fill!(copy(bstates_sparse2), sqrt(ΔE_sparse2 / (M_sparse-1)))
    bathcoupling = [bcoupling_sparse1; bcoupling_window; bcoupling_sparse2]

    return WindowedTrapezoidalRule(bathstates, bathcoupling)
end