
module SchNetPackModels

# Attention
# in order to make spk work in julia, you need to change the following line after installation
# in src/schnetpack/nn/cfconv.py change line 78 to "nbh = neighbors.reshape(1,nbh_size[1]*nbh_size[2],1)"

using DocStringExtensions
using UnPack
using ..PyCall
using ....NonadiabaticMolecularDynamics: AbstractCell, PeriodicCell, Atoms
using ..Models
using PeriodicTable
using Unitful, UnitfulAtomic

export SchNetPackModel
export FrictionSchNetPackModel

const torch = PyNULL() # Julia torch module
const ase = PyNULL()
const spk = PyNULL()

function __init__()
    copy!(torch, pyimport_conda("torch", "torch"))
    copy!(spk, pyimport_conda("schnetpack", "schnetpack", "conda-forge"))
    copy!(ase, pyimport_conda("ase", "ase=3.21.1", "conda-forge"))
end

include("ML_descriptor.jl")

"""
$(TYPEDEF)

Internal type for storing parameters, inputs and the model from SchNetPack.
"""
struct SchNetPackInterface{C<:AbstractCell,A<:Atoms}
    cell::C
    atoms::A
    input::Dict{String,PyObject}
    args::PyObject
    units::Dict{String,Unitful.Units}
    model::PyObject
end

function SchNetPackInterface(path, atoms, cell)
    model = torch_load(path)
    args = get_model_args(path)
    units = get_units(path, args)

    args.atomic_charges = atoms.numbers
    args.pbc = cell.periodicity
    
    input = Dict{String, PyObject}()
    SchNetPackInterface(cell, atoms, input, args, units, model)
end

"Load model using torch from file."
function torch_load(path::String)
    device = "cpu"
    torch.load(joinpath(path, "best_model"), map_location=torch.device(device)).to(device)
end

"Read the arguments for the model from `args.json`."
function get_model_args(path::String)
    args = spk.utils.read_from_json(joinpath(path,"args.json"))
    if args.cuda == true && torch.cuda.is_available() == true
        args.device = "cuda"
    else 
        args.device = "cpu"
    end

    if args.force === nothing && args.negative_dr == false
        args.sign = 1 # Positive sign if trained on gradients.
    else
       args.sign = -1 # Negative if trained on forces.
    end

    return args
end

"Extract the units from the model args and convert to `Unitful.Units`."
function get_units(path::String, args::PyObject)::Dict{String,Unitful.Units}
    dataset = spk.data.AtomsData(string(joinpath(path, args.datapath)),collect_triples=args=="wacsf")
    if in("units",keys(dataset.get_metadata()))
        units = dataset.get_metadata("units")
        map!(u->uparse(u), values(units))
        map!(u-> u isa Unitful.Units ? u : unit(u), values(units))
    else
        println("Attention! No units found in the data set. Atomic units are used.")
        units = Dict()
        units["energy"] = u"hartree"
        units["forces"] = u"hartree/bohr"
        #TODO check EFT units
        units["EFT"] = u"hartree/bohr^2"
        units["R"] = u"bohr"
    end
    return units
end

"""
$(TYPEDEF)

Model for adiabatic dynamics with a SchNetPack trained model.
"""
struct SchNetPackModel{I} <: Models.AdiabaticModel
    pes_interface::I
end

function SchNetPackModel(path::String, cell::PeriodicCell, atoms::Atoms)
    pes_interface = SchNetPackInterface(path, atoms, cell)
    SchNetPackModel(pes_interface)
end

"""
$(TYPEDEF)

Model for MDEF with a SchNetPackFriction model.
"""
struct FrictionSchNetPackModel{M,I2} <: AdiabaticFrictionModel
    pes_model::M
    friction_interface::I2
    F_idx::Vector{Int}
end

function FrictionSchNetPackModel(pes_model, friction_path, cell, atoms, F_idx)
    friction_interface = SchNetPackInterface(friction_path, atoms, cell)
    FrictionSchNetPackModel(pes_model, friction_interface, F_idx)
end

function Models.potential!(model::SchNetPackModel, V::AbstractVector, R::AbstractMatrix)
    @unpack input, cell, atoms, args, units = model.pes_interface
    update_schnet_input!(input, cell, atoms, R, args, units)
    V .= austrip.(model.pes_interface.model(input)["energy"].detach().numpy()[1] .* units["energy"])
end

function Models.derivative!(model::SchNetPackModel, D::AbstractMatrix, R::AbstractMatrix)
    @unpack input, cell, atoms, args, units = model.pes_interface
    update_schnet_input!(input, cell, atoms, R, args, units)
    sign = args.sign # spk already multiplies with -1 if forces are trained
    D .= austrip.(sign*model.pes_interface.model(input)["forces"].detach().numpy()[1,:,:]' .* units["forces"])
end

Models.potential!(model::FrictionSchNetPackModel, V, R) = Models.potential!(model.pes_model, V, R)
Models.derivative!(model::FrictionSchNetPackModel, D, R) = Models.derivative!(model.pes_model, D, R)

function Models.friction!(model::FrictionSchNetPackModel, F::AbstractMatrix, R::AbstractMatrix)
    @unpack input, cell, atoms, args, units = model.friction_interface
    update_schnet_input_friction!(input, cell, atoms, R, args, units, model.F_idx)
    friction = model.friction_interface.model(input)["friction_tensor"].detach().numpy()
    indices = model.F_idx
    @views for (m,j) in enumerate(indices)
        for (n,i) in enumerate(indices)
            F[3(i-1)+1:3i, 3(j-1)+1:3j] .= austrip.(friction[1,3(n-1)+1:3n, 3(m-1)+1:3m] .* units["EFT"])
        end
    end

    masses = model.friction_interface.atoms.masses
    for j in indices
        for i in indices
            F[3(i-1)+1:3i, 3(j-1)+1:3j] .*= sqrt(masses[i] * masses[j])
        end
    end

    return F
end

end # module
