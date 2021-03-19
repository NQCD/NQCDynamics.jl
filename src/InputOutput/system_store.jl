using FileIO

export SystemStore
export save
export load

"""
$(TYPEDEF)

Container for storing system information.

This can be used to store simulation results or configurations between
sampling and running dynamics.
"""
struct SystemStore{T,A<:Atoms,C<:AbstractCell}
    data::T
    atoms::A
    cell::C
end

save(f, system::SystemStore) = FileIO.save(f, "system", system)
load(f) = FileIO.load(f, "system")
