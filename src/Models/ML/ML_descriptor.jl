using Unitful
using PyCall
using ....Atoms

const torch = PyNULL() # Julia torch module

export update_schnet_input!

# pushfirst!(PyVector(pyimport("sys")."path"),"")

function update_schnet_input!(schnet_inputs::Dict, p::AtomicParameters, R::AbstractVector, model_args::PyObject)

    cell = ustrip.(u"Å", p.cell.vectors)
    positions = reshape(R, (p.n_atoms, 3))
    # We might want to get Julia versions of these, also we should remove the if statement.
    if model_args.environment_provider == "simple"
        nbh_idx, offsets = py"get_simple_environment"(p.n_atoms,p.atom_types,positions,model_args.pbc,cell,model_args)
    elseif model_args.environment_provider == "ase"
        nbh_idx,offsets = py"get_ase_environment"(p.n_atoms,p.atom_types,positions,model_args.pbc,cell,model_args)
    elseif model_args.environment_provider == "torch"
        nbh_idx,offsets = py"get_torch_environment"(p.n_atoms,model_args.atomic_charges,positions,model_args.pbc,cell,model_args)
    end
    # Some of these will not change and do not need to be updated every time.
    schnet_inputs["_atomic_numbers"] =  torch.LongTensor(model_args.atomic_charges).unsqueeze(0).to(model_args.device)
    if model_args.contributions == true
        #manual for 2 atoms
        #TODO update if we ever use larger molecules 
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
    schnet_inputs
end


function __init__()
    copy!(torch, pyimport("torch"))
    py"""
    import schnetpack as spk
    from schnetpack import Properties
    import torch
    import numpy as np
    import ase 
    import math

    #COPYRIGHT SPK: adapted to work without ase atoms

    def get_simple_environment(n_atoms, species_, positions_,pbc_array,cell_,model_args):
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

    def get_ase_environment(n_atoms, species_, positions_,pbc_array,cell_,model_args):
        cutoff = model_args.cutoff
        atoms=ase.Atoms(species_,np.array(positions_))
        atoms.set_cell(cell_)
        atoms.set_pbc(pbc_array)
        idx_i, idx_j, idx_S = ase.neighborlist.neighbor_list(
            "ijS", atoms, cutoff, self_interaction=False
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

    def get_torch_environment(n_atoms, species_, positions_,pbc_array,cell_,model_args):
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

    def compute_shifts(cell, pbc, cutoff):

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

    def neighbor_pairs(padding_mask, coordinates, cell, shifts, cutoff):
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