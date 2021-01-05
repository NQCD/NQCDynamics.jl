module ML

using PyCall
using PeriodicTable
using ....Atoms
using ..Models
using UnitfulAtomic
export SchNetPackModel
export initialize_MLmodel
export torch_load
export getmetadata
export get_param_model
export check_atoms
np = pyimport("numpy")
ase = pyimport("ase")

include("ML_descriptor.jl")

struct SchNetPackModel <: Models.Model

    n_states::UInt
    get_V0::Function
    get_D0::Function

    function SchNetPackModel(path::String, atoms::AtomicParameters)
        model, model_args, ML_units = initialize_MLmodel(path, atoms)
        input = Dict{String, PyObject}()
        energy_unit = austrip(1*uparse(ML_units["energy"]))
        R_unit = 1/austrip(1*uparse(ML_units["R"]))
        force_unit = energy_unit*R_unit
        #TODO Check EFT unit
        EFT_unit = energy_unit/R_unit/R_unit

        #Comment/TODO: the schnet update only needs to be done once for energy and is not needed for the forces.
        #for eft the schnet update needs to be done as less atoms are needed
        function get_V0(R::AbstractVector)
            update_schnet_input!(input, atoms, (R*R_unit), model_args)
            model(input)["energy"].detach().numpy()[1,1]*energy_unit
        end

        function get_D0(R::AbstractVector)::Vector
            update_schnet_input!(input, atoms, (R/(austrip(1*uparse(ML_units["R"])))), model_args)
            -model(input)["forces"].detach().numpy()[1,:,:]'[:]*force_unit
        end
        
        new(1, get_V0, get_D0)
    end
end

function initialize_MLmodel(path::String, atoms::AtomicParameters)
    model = torch_load(path)
    model_args = get_param_model(path)
    #get metadata
    ML_units = get_metadata(path, model_args, atoms.atom_types)
    force_mask = false
    #global for now
    model_args.atomic_charges = [elements[type].number for type in atoms.atom_types]
    model_args.pbc = atoms.cell.periodicity
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

# This function is not type stable. force_mask can be either Bool or Vector
# This should be removed
function get_metadata(path::String, model_args::PyObject, atoms::Vector{Symbol})
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