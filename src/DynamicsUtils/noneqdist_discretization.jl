using Random,Distributions,DelimitedFiles,DataInterpolations,FastGaussQuadrature,Trapz # removed Plots from the using list

"""
Discretization script generates (in this default case) 1000 occupation vectors denoting which states in the discretized
density of states are occupied. 

The variable `binary_vectors` / `splits` contains the 1000 discretized versions of the non-equilibrium distributions i.e. our 1000 occupation vectors.
The NAH density of states for the bath has been taken to have 250 discretized states, hence each occupation vector is 250 elements long.

The "average" of these 1000 occupation vectors return the original non-equilibrium distribution. 
"""

# -------------------------------------- Necesary Functions -------------------------------------- #
function generate_DOS(File::String,n) #Simply builds a spline of a DOS from a file and scales from eV^-1atom^-1 to eV^-1nm^-1
    TotalDOS::Matrix{Float64}=readdlm(File)
    return DataInterpolations.LinearInterpolation(TotalDOS[:,2].*n,TotalDOS[:,1],extrapolate=true)
end

function grid_builder(l,Espan) #Builds a grid similar to that used in News Andersen Hamiltonian - redundant
    gh,weights = gausshermite(l)
    return ((gh .- minimum(gh)) ./ (maximum(gh)/(Espan/2))) .- (Espan/2) 
end

function discretization(dis, n::Int, dos,grid;mean_tol=0.005,particle_tol=0.02) 
    #Converts a vector x into n binary vectors such that their mean equals x and idnvidually maintain particle conservation given by w
    @assert length(dis) == length(dos) "x and w must have the same length"
    
    m = length(dis)
    binary_vectors = zeros(Int, m, n) #Generates matrix for our binary distribuutions

    for i in 1:m
        for j in 1:n
            binary_vectors[i, j] = rand(Bernoulli(dis[i])) #Fills with random numbers of 0 or 1
        end
    end

    particle_number = trapz(grid,dos.*dis) #Calculates total number of particles in original distribution using uenven trapezium rule
    tolerance_particle = particle_number * particle_tol  # Tolerance for particle conservation

    for j in 1:n
        current_particles = trapz(grid, binary_vectors[:, j] .* dos)  # Particles in current binary vector
        diff = particle_number - current_particles                    # Difference in particle count
        
        while abs(diff) > tolerance_particle  # Continue adjusting until within tolerance
            if diff > 0
                # We need to increase particles -> change 0's to 1's
                for i in 1:m
                    if abs(diff) <= tolerance_particle
                        break
                    elseif binary_vectors[i, j] == 0 && dos[i] > 0
                        # Change 0 to 1, but ensure that we don't violate the mean condition
                        current_mean = sum(binary_vectors[i, :]) / n
                        if current_mean + (1 / n) - dis[i] < mean_tol  # Add tolerance for mean checking
                            binary_vectors[i, j] = 1
                            current_particles = trapz(grid, binary_vectors[:, j] .* dos)
                            diff = particle_number - current_particles
                        end
                    end
                end
            elseif diff < 0
                # We need to decrease particles -> change 1's to 0's
                for i in m:-1:1
                    if abs(diff) <= tolerance_particle
                        break
                    elseif binary_vectors[i, j] == 1 && dos[i] > 0
                        # Change 1 to 0, but ensure that we don't violate the mean condition
                        current_mean = sum(binary_vectors[i, :]) / n
                        if current_mean - (1 / n) - dis[i] > -mean_tol  # Add tolerance for mean checking
                            binary_vectors[i, j] = 0
                            current_particles = trapz(grid, binary_vectors[:, j] .* dos)
                            diff = particle_number - current_particles
                        end
                    end
                end
            end
        end
    end

    cps=zeros(n)
    for j in 1:n
        cps[j] = trapz(grid,binary_vectors[:, j] .* dos)
    end
    pushfirst!(cps,particle_number)

    return binary_vectors,cps
end

function restore_ftot(splits,n) #Re-averages the binary dsitributions to check they rebuild true distribution
    restore=zeros(length(splits[:,1]))
    for i in eachindex(splits[:,1])
        restore[i]=sum(splits[i,:])
    end
    return restore./n
end
# ------------------------------------------------------------------------------------------------ #


# # ------------------------- Generation of discretized Neq dist. over DOS ------------------------- #
# dis_from_file = readdlm("Neq_distribution.csv") # Read distribution from file
# energy_grid = dis_from_file[:,1] # Seperate into energy grid 
# distribution = dis_from_file[:,2] # and total distribution
# dis_spl = LinearInterpolation(distribution,energy_grid) # Spline the distribution to project onto new energy grid
# DOS_spl = generate_DOS("Cu_DOS.dat",85) # Build DOS spline from file

# NAH_energygrid = grid_builder(250,10.0) # Build new energy grid
# DOS = DOS_spl(NAH_energygrid) # Recast DOS and distribution onto new grid
# dis = dis_spl(NAH_energygrid)

# vecs=1000 # Number of binary vectors requested
# splits,parts = discretization(dis,vecs,DOS,NAH_energygrid,mean_tol=0.01,particle_tol=0.001) #Builds the vectors and returns the vectors (splits) and the relative particle distribution (parts)
# #The first value of parts is the correct value if you would like to plot a h-line or something to see the distribution

# # ----------------------------- Regenerate original Neq distribution ----------------------------- #
# restored_ftot = restore_ftot(splits,vecs) # Rebuilds original distribution