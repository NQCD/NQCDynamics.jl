using NonadiabaticMolecularDynamics

n_DoF = 3
n_atoms = 20

# Can construct using matrix with each atom as columns 
R = Positions(rand(n_DoF, n_atoms))
P = Momenta(rand(n_DoF, n_atoms))

# Alternatively, construct with array of arrays
R = Positions([rand(n_DoF) for i=1:n_atoms])
P = Momenta([rand(n_DoF) for i=1:n_atoms])

R[:, 5] # Access atom 5
R[5] # Access atom 5
P[:, 5] # Access atom 5
P[5] # Access atom 5

# Standard array operations should work
R + P

using RecursiveArrayTools
typeof(ArrayPartition(R, P))
a = Phasespace(R, P)
typeof(a.x.x[1])
