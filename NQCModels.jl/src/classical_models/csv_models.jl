using BSplineKit

#==
Since CSV files come in many shapes and sizes, leave importing them into a matrix to the user. 

CSV model needs an input matrix of shape (n,2), where the first column of the matrix denotes the position, the second half the potential value in Hartree. 
The model assumes that the positions are evenly spaced, forming a uniform grid, for faster fitting. 
==#

```
This package contains a type of NQCModel which takes a Matrix of positions and potential values as an input to fit a spline to. 
Currently, only a one-dimensional adiabatic model is implemented: CSVModel_1D Importing the table of positions and potential values is 
left to the user since CSV files come in too many different flavours.
```
Parameters.@with_kw struct CSVModel_1D <: ClassicalModel
    # Input matrix
    potential_matrix::Matrix{Float64}
    # build interpolation and derivative object once, then evaluate with functions
    potential_function=interpolate(potential_matrix[:,1],potential_matrix[:,2],BSplineOrder(3))
    derivative_function=Derivative(1)*potential_function
end



NQCModels.ndofs(model::CSVModel_1D)=1

function NQCModels.potential(model::CSVModel_1D, R::AbstractMatrix)
    return(model.potential_function(R[1,1]))
end

function NQCModels.derivative!(model::CSVModel_1D, D::AbstractMatrix, R::AbstractMatrix)
    D.=model.derivative_function(R[1,1])
end

CSVModel_1D(x)=CSVModel_1D(potential_matrix=x) # Definition shortcut. 

