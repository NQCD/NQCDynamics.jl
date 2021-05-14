
function update_schnet_input!(schnet_inputs::Dict, periodic_cell::PeriodicCell, atoms::Atoms, R::AbstractMatrix, model_args::PyObject, ML_units::Dict)
    cell = ustrip.(auconvert.(u"Å", periodic_cell.vectors))
    #cutoff radius should always be in the correct units - this should not be changed
    positions = ustrip.(auconvert.(ML_units["R"], R'))
    molecule = ase.atoms.Atoms(atoms.types,positions)
    molecule.set_cell(cell)
    molecule.set_pbc(model_args.pbc)
    # We might want to get Julia versions of these, also we should remove the if statement.
    if model_args.environment_provider == "simple"
        nbh_idx, offsets = spk.environment.SimpleEnvironmentProvider().get_environment(molecule)
        #py"get_simple_environment"(length(atoms), atoms.types, positions, model_args.pbc, cell, model_args)
    elseif model_args.environment_provider == "ase"
        nbh_idx,offsets = spk.environment.AseEnvironmentProvider(model_args.cutoff).get_environment(molecule)
        #py"get_ase_environment"(length(atoms), atoms.types, positions, model_args.pbc, cell, model_args)
    elseif model_args.environment_provider == "torch"
        nbh_idx,offsets = spk.environment.TorchEnvironmentProvider(model_args.cutoff,model_args.device).get_environment(molecule)
        #py"get_torch_environment"(length(atoms), model_args.atomic_charges, positions, model_args.pbc, cell, model_args)
    end
    schnet_inputs["_atomic_numbers"] =  torch.LongTensor(molecule.get_atomic_numbers()).unsqueeze(0).to(model_args.device)
    if model_args.contributions == true
        schnet_inputs["_atom_mask"] = torch.zeros_like(schnet_inputs["_atomic_numbers"]).float()
        schnet_inputs["_atom_mask"][0,2].fill_(1.0)
        schnet_inputs["_atom_mask"][0,1].fill_(1.0)
        schnet_inputs["_atom_mask"].unsqueeze(0).to(model_args.device)
    else
        schnet_inputs["_atom_mask"] = torch.ones_like(schnet_inputs["_atomic_numbers"]).float().unsqueeze(0).to(model_args.device)
    end
    schnet_inputs["_positions"] = torch.FloatTensor(positions).unsqueeze(0).to(model_args.device)
    mask = torch.FloatTensor(nbh_idx) >= 0
    schnet_inputs["_neighbor_mask"] = mask.float().unsqueeze(0).to(model_args.device)
    schnet_inputs["_neighbors"] = (torch.LongTensor(nbh_idx) * mask.long()).unsqueeze(0).to(model_args.device)
    schnet_inputs["_cell"] = torch.FloatTensor(cell).unsqueeze(0).to(model_args.device)
    schnet_inputs["_cell_offset"] = torch.FloatTensor(offsets).unsqueeze(0).to(model_args.device).contiguous()
end

function update_schnet_input_friction!(schnet_inputs::Dict, periodic_cell::PeriodicCell, atoms::Atoms, R::AbstractMatrix, model_args::PyObject, ML_units::Dict,friction_indices::Vector{Int64})
    
    update_schnet_input!(schnet_inputs, periodic_cell, atoms, R, model_args, ML_units)
    schnet_inputs["friction_indices"] = torch.zeros((1,length(friction_indices)))
    it=0
    for i in friction_indices
        it+=1
        schnet_inputs["friction_indices"][0,i+1]=Int8(friction_indices[it])
    end
    schnet_inputs["_atom_mask"] = schnet_inputs["_atom_mask"].reshape(1,-1)
    schnet_inputs["friction_indices"]=schnet_inputs["friction_indices"].int()
    schnet_inputs["_idx"]=torch.Tensor(1,1).int()
    #torch_friction_indices = schnet_inputs["friction_indices"][0].int()
    #friction_input = torch.index_select(schnet_inputs["_positions"],1,torch_friction_indices)
    #friction_input.requires_grad = true
    #julia starts to count from 1
    #it=0
    #for i in friction_indices
    #    it+=1
    #    schnet_inputs["_positions"][0,i+1] = friction_input[0,it]
    #end
    #setindex does not work with pyobjects
    #schnet_inputs["_positions"][:,friction_indices,:] .= friction_input
    #schnet_inputs["friction_positions"] = friction_input

end
