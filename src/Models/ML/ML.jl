module ML

using PyCall
using PeriodicTable
using ....Atoms
using ..Models

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
        model, model_args, force_mask = initialize_MLmodel(path, atoms)
        input = Dict{String, PyObject}()
        
        function get_V0(R::AbstractVector)
            update_schnet_input!(input, atoms, R, model_args)
            model(input)["energy"].detach().numpy()[1,1]
            print(model(input))
        end

        function get_V0_atomic(R::AbstractVector)
            update_schnet_input!(input, atoms, R, model_args)
        end
        
        function get_D0(R::AbstractVector)::Vector
            update_schnet_input!(input, atoms, R, model_args)
            -model(input)["forces"].detach().numpy()[1,:,:]'[:]
        end
        
        new(1, get_V0, get_D0)
    end
end

function initialize_MLmodel(path::String, atoms::AtomicParameters)
    model = torch_load(path)
    model_args = get_param_model(path)
    #get metadata
    #force_mask = get_metadata(path, model_args, atoms.atom_types)
    force_mask = false
    #global for now
    model_args.atomic_charges = [elements[type].number for type in atoms.atom_types]
    model_args.pbc = atoms.cell.periodicity
    return model, model_args, force_mask
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

    if in("force_mask",keys(dataset.get_metadata()))
        index = length(dataset.get_metadata("force_mask_index")) - length(atoms) 
        if index <= 0
            return dataset.get_metadata("force_mask_index")
        else 
            return dataset.get_metadata("force_mask_index")[index:length(dataset.get_metadata("force_mask_index"))]
        end
    else
        return force_mask
    end
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