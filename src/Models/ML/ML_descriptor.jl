module ML_descriptor

using PyCall
#import the necessary functions
np = pyimport("numpy")
ase = pyimport("ase")
spk = pyimport("schnetpack")
using Unitful
using UnitfulAtomic    
export julia2schnet
export schnet2julia
pushfirst!(PyVector(pyimport("sys")."path"),"")

#abstract type DynamicalVariables end

#struct Descriptor_Variables{T<:AbstractFloat} <: DynamicalVariables
#    schnet_inputs::Dict
#end
function pyinit()
    py"""
    import schnetpack as spk
    from schnetpack import Properties
    import torch
    import numpy as np
    import ase 
    from ase import neighborlist


    #COPYRIGHT SPK: adapted to work without ase atoms



    def get_simple_environment(n_atoms, species_, positions_,pbc_array,cell_,model_args):
            import numpy as np
            if n_atoms == 1:
                neighborhood_idx = -np.ones((1, 1), dtype=np.float32)
                offsets = np.zeros((n_atoms, 1, 3), dtype=np.float32)
            else:
                neighborhood_idx = np.tile(
                    np.arange(n_atoms, dtype=np.float32)[np.newaxis], (n_atoms, 1)
                )

                neighborhood_idx = neighborhood_idx[
                    ~np.eye(n_atoms, dtype=np.bool)
                ].reshape(n_atoms, n_atoms - 1)

                offsets = np.zeros(
                    (neighborhood_idx.shape[0], neighborhood_idx.shape[1], 3),
                    dtype=np.float32,
                )
            return neighborhood_idx, offsets
    """

    py"""
    def get_ase_environment(n_atoms, species_, positions_,pbc_array,cell_,model_args):
            import numpy as np
            from ase.neighborlist import neighbor_list
            cutoff = model_args.cutoff

            idx_i, idx_j, idx_S = neighbor_list(
                "ijS", ase.Atoms(species_,positions_), cutoff, self_interaction=False
            )
            if idx_i.shape[0] > 0:
                uidx, n_nbh = np.unique(idx_i, return_counts=True)
                n_max_nbh = np.max(n_nbh)

                n_nbh = np.tile(n_nbh[:, np.newaxis], (1, n_max_nbh))
                nbh_range = np.tile(
                    np.arange(n_max_nbh, dtype=np.int)[np.newaxis], (n_nbh.shape[0], 1)
                )

                mask = np.zeros((n_atoms, np.max(n_max_nbh)), dtype=np.bool)
                mask[uidx, :] = nbh_range < n_nbh
                neighborhood_idx = -np.ones((n_atoms, np.max(n_max_nbh)), dtype=np.float32)
                neighborhood_idx[mask] = idx_j

                offset = np.zeros((n_atoms, np.max(n_max_nbh), 3), dtype=np.float32)
                offset[mask] = idx_S
            else:
                neighborhood_idx = -np.ones((n_atoms, 1), dtype=np.float32)
                offset = np.zeros((n_atoms, 1, 3), dtype=np.float32)

            return neighborhood_idx, offset
    """
    py"""
    def get_torch_environment(n_atoms, species_, positions_,pbc_array,cell_,model_args):
            import torch 
            import numpy as np
            device = str(model_args.device)
            species = torch.FloatTensor(species_).to(device)
            coordinates = torch.FloatTensor(positions_).to(device)
            pbc = torch.from_numpy(pbc_array.astype("uint8")).to(device)
            cutoff = model_args.cutoff 
            
            if not cell_.any():
                cell = torch.eye(3, dtype=species.dtype).to(device)
            else:
                cell = torch.Tensor(cell_).to(device)
            shifts = compute_shifts(cell=cell,pbc=pbc,cutoff=cutoff)
            print(shifts)


            # The returned indices are only one directional
            idx_i, idx_j, idx_S = neighbor_pairs(
                species == -1, coordinates, cell, shifts, cutoff
            )


            idx_i = idx_i.cpu().detach().numpy()
            idx_j = idx_j.cpu().detach().numpy()
            idx_S = idx_S.cpu().detach().numpy()

            # Create bidirectional id arrays, similar to what the ASE neighbor_list returns
            bi_idx_i = np.hstack((idx_i, idx_j))
            bi_idx_j = np.hstack((idx_j, idx_i))
            bi_idx_S = np.vstack((-idx_S, idx_S))

            if bi_idx_i.shape[0] > 0:
                uidx, n_nbh = np.unique(bi_idx_i, return_counts=True)
                n_max_nbh = np.max(n_nbh)

                n_nbh = np.tile(n_nbh[:, np.newaxis], (1, n_max_nbh))
                nbh_range = np.tile(
                    np.arange(n_max_nbh, dtype=np.int)[np.newaxis], (n_nbh.shape[0], 1)
                )

                mask = np.zeros((n_atoms, np.max(n_max_nbh)), dtype=np.bool)
                mask[uidx, :] = nbh_range < n_nbh
                neighborhood_idx = -np.ones((n_atoms, np.max(n_max_nbh)), dtype=np.float32)
                offset = np.zeros((n_atoms, np.max(n_max_nbh), 3), dtype=np.float32)

                # Assign neighbors and offsets according to the indices in bi_idx_i, since in contrast
                # to the ASE provider the bidirectional arrays are no longer sorted.
                # TODO: There might be a more efficient way of doing this than a loop
                for idx in range(n_atoms):
                    neighborhood_idx[idx, mask[idx]] = bi_idx_j[bi_idx_i == idx]
                    offset[idx, mask[idx]] = bi_idx_S[bi_idx_i == idx]

            else:
                neighborhood_idx = -np.ones((n_atoms, 1), dtype=np.float32)
                offset = np.zeros((n_atoms, 1, 3), dtype=np.float32)

            return neighborhood_idx, offset
    """
    py"""
    def compute_shifts(cell, pbc, cutoff):
        import math
        import torch

        # type: (Tensor, Tensor, float) -> Tensor
        reciprocal_cell = cell.inverse().t()
        inv_distances = reciprocal_cell.norm(2, -1)
        num_repeats = torch.ceil(cutoff * inv_distances).to(torch.long)
        num_repeats = torch.where(pbc, num_repeats, torch.zeros_like(num_repeats))

        r1 = torch.arange(1, num_repeats[0] + 1, device=cell.device)
        r2 = torch.arange(1, num_repeats[1] + 1, device=cell.device)
        r3 = torch.arange(1, num_repeats[2] + 1, device=cell.device)
        o = torch.zeros(1, dtype=torch.long, device=cell.device)

        return torch.cat(
            [
                torch.cartesian_prod(r1, r2, r3),
                torch.cartesian_prod(r1, r2, o),
                torch.cartesian_prod(r1, r2, -r3),
                torch.cartesian_prod(r1, o, r3),
                torch.cartesian_prod(r1, o, o),
                torch.cartesian_prod(r1, o, -r3),
                torch.cartesian_prod(r1, -r2, r3),
                torch.cartesian_prod(r1, -r2, o),
                torch.cartesian_prod(r1, -r2, -r3),
                torch.cartesian_prod(o, r2, r3),
                torch.cartesian_prod(o, r2, o),
                torch.cartesian_prod(o, r2, -r3),
                torch.cartesian_prod(o, o, r3),
            ]
        )
    """
    py"""
    def neighbor_pairs(padding_mask, coordinates, cell, shifts, cutoff):
        import math
        import torch
        # type: (Tensor, Tensor, Tensor, Tensor, float) -> Tuple[Tensor, Tensor, Tensor, Tensor]

        coordinates = coordinates.detach()
        cell = cell.detach()
        num_atoms = padding_mask.shape[0]
        all_atoms = torch.arange(num_atoms, device=cell.device)

        # Step 2: center cell
        p12_center = torch.triu_indices(num_atoms, num_atoms, 1, device=cell.device)
        shifts_center = shifts.new_zeros((p12_center.shape[1], 3))

        # Step 3: cells with shifts
        # shape convention (shift index, molecule index, atom index, 3)
        num_shifts = shifts.shape[0]
        all_shifts = torch.arange(num_shifts, device=cell.device)
        prod = torch.cartesian_prod(all_shifts, all_atoms, all_atoms).t()
        shift_index = prod[0]
        p12 = prod[1:]
        shifts_outside = shifts.index_select(0, shift_index)

        # Step 4: combine results for all cells
        shifts_all = torch.cat([shifts_center, shifts_outside])
        p12_all = torch.cat([p12_center, p12], dim=1)
        shift_values = shifts_all.to(cell.dtype) @ cell

        # step 5, compute distances, and find all pairs within cutoff
        selected_coordinates = coordinates.index_select(0, p12_all.view(-1)).view(2, -1, 3)
        distances = (
            selected_coordinates[0, ...] - selected_coordinates[1, ...] + shift_values
        ).norm(2, -1)
        padding_mask = padding_mask.index_select(0, p12_all.view(-1)).view(2, -1).any(0)
        distances.masked_fill_(padding_mask, math.inf)
        in_cutoff = (distances < cutoff).nonzero()

        pair_index = in_cutoff.squeeze()
        atom_index12 = p12_all[:, pair_index]
        shifts = shifts_all.index_select(0, pair_index)
        atom_index1 = atom_index12[0, ...]
        atom_index2 = atom_index12[1, ...]

        return atom_index1, atom_index2, shifts
    """
end
#abstract type DynamicalVariables end

#s#truct schnet_inputs_julia{T<:AbstractFloat} <: DynamicalVariables
# #   atomic_charges::PyArray
#end



function julia2schnet(p::Any,R::Array,model_args::PyObject)
    torch = pyimport("torch")
    spk=pyimport("schnetpack")

    cell = ustrip.(u"Å", p.cell.vectors)
    positions = positions_spk(R,p.n_atoms)
    pyinit()
    if model_args.environment_provider == "simple"
        #nbh_idx, offsets = simple_env(p)
        nbh_idx, offsets = py"get_simple_environment"(p.n_atoms,p.atom_types,positions,model_args.pbc,cell,model_args)
    elseif model_args.environment_provider == "ase"
        nbh_idx,offsets = py"get_ase_environment"(p.n_atoms,p.atom_types,positions,model_args.pbc,cell,model_args)
    elseif model_args.environment_provider == "torch"
        nbh_idx,offsets = py"get_torch_environment"(p.n_atoms,model_args.atomic_charges,positions,model_args.pbc,cell,model_args)
    end
    
    schnet_inputs = Dict()
    schnet_inputs[spk.Properties.Z] =  torch.LongTensor(model_args.atomic_charges)
    schnet_inputs[spk.Properties.atom_mask] = torch.ones_like(schnet_inputs[spk.Properties.Z]).float()
    schnet_inputs[spk.Properties.R] = torch.FloatTensor(positions)
    mask = torch.FloatTensor(nbh_idx) >= 0
    schnet_inputs[spk.Properties.neighbor_mask] = mask.float()
    schnet_inputs[spk.Properties.neighbors] = torch.LongTensor(nbh_idx) * mask.long()
    #get cells
    schnet_inputs[spk.Properties.cell] = torch.FloatTensor(cell)
    schnet_inputs[spk.Properties.cell_offset]= torch.FloatTensor(offsets)
    schnet_inputs_julia=Dict()
    for (k,v) in schnet_inputs
        schnet_inputs_julia[k] = v.unsqueeze(0).to(model_args.device)
    end
    return schnet_inputs_julia
end

function positions_spk(R,n_atoms)
    return reshape(R,(n_atoms,3))
end

function schnet2julia()

end
    

end