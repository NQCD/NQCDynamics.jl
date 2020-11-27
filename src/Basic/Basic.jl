"""
Data containers and functions for preparing the system. 
"""
module Basic

using RecursiveArrayTools
using PeriodicTable
using Unitful
using UnitfulAtomic
using PyCall

include("cell.jl")
include("system.jl")
include("phasespace.jl")

include("io.jl")

end # module


