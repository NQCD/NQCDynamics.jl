"""
Data containers and functions for preparing the system. 
"""
module Basic

using RecursiveArrayTools
using PeriodicTable
using Unitful
using UnitfulAtomic

include("cell.jl")
include("system.jl")
include("phasespace.jl")

end # module


