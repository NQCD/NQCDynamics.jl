module ML

using PyCall
using PeriodicTable
using ....NonadiabaticMolecularDynamics
using ..Models
using UnitfulAtomic
using Unitful

export SchNetPackModel

np = pyimport("numpy")
ase = pyimport("ase")

include("ML_descriptor.jl")

struct SchNetPackModel <: AdiabaticModel
    cell::PeriodicCell
    atoms::Atoms
    input::Dict{String, PyObject}
    model_args::PyObject
    ML_units::Dict{Any, Any}
    model::PyObject
end

function SchNetPackModel(path::String, cell::PeriodicCell, atoms::Atoms)
    model, model_args, ML_units = initialize_MLmodel(path, cell, atoms)
    input = Dict{String, PyObject}()
    SchNetPackModel(cell, atoms, input, model_args, ML_units, model)
end

function Models.potential!(model::SchNetPackModel, V::AbstractVector, R::AbstractMatrix)
    update_schnet_input!(model.input, model.cell, model.atoms, R, model.model_args, model.ML_units)
    V .= austrip.(model.model(model.input)["energy"].detach().numpy()[1] .* uparse(model.ML_units["energy"]))
end

function Models.derivative!(model::SchNetPackModel, D::AbstractMatrix, R::AbstractMatrix)
    update_schnet_input!(model.input, model.cell, model.atoms, R, model.model_args, model.ML_units)
    D .= austrip.(-model.model(model.input)["forces"].detach().numpy()[1,:,:]' .* uparse(model.ML_units["forces"]))
end

function initialize_MLmodel(path::String, cell::PeriodicCell, atoms::Atoms)
    model = torch_load(path)
    model_args = get_param_model(path)
    #get metadata
    ML_units = get_metadata(path, model_args)
    force_mask = false
    model_args.atomic_charges = [elements[type].number for type in atoms.types]
    model_args.pbc = cell.periodicity
    return model, model_args, ML_units
end

function torch_load(path::String)
    torch = pyimport("torch")
    device = "cpu"
    torch.load(joinpath(path, "best_model"), map_location=torch.device(device)).to(device)
end

function get_param_model(path::String)
    spk = pyimport("schnetpack")
    model_args = spk.utils.read_from_json(joinpath(path,"args.json"))
    if model_args.cuda == true && torch.cuda.is_available() == true
        model_args.device = "cuda"
    else 
        model_args.device = "cpu"
    end
    #environment_provider = spk.utils.script_utils.settings.get_environment_provider(model_args,device=device)
    return model_args #, environment_provider
end

function get_metadata(path::String, model_args::PyObject)
    spk = pyimport("schnetpack")
    dataset = spk.data.AtomsData(string(joinpath(path,model_args.datapath)),collect_triples=model_args=="wacsf")
    if in("units",keys(dataset.get_metadata()))
        ML_units = dataset.get_metadata("units")
    else
        println("Attention! No units found in the data set. Atomic units are used.")
        ML_units = Dict()
        ML_units["energy"] = "Eₕ"
        ML_units["forces"] = "Eₕ/a₀"
        #TODO check EFT units
        ML_units["EFT"] = "Eₕ/a₀/a₀"
        ML_units["R"] = "a₀"
    end
    return ML_units
end

function get_properties(model_path, atoms, args...)
    #get energy and forces and other relevant properties
    #only for energy and properties with simple spk
    atoms.set_calculator(calculator)
    energy = atoms.get_total_energy()
    if force_mask == false
        forces = atoms.get_forces()
    else
        forces = atoms.get_forces() * force_mask
    end

    return (energy,forces)
end

end # module