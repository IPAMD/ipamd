"""
plugin for calculating the radial distribution function g(r)
"""
import math
import numpy as np
from numba import cuda
from ipamd.public.models.data import Vector
from ipamd.public.utils.parser import protein_range_split

configure = {
    "schema": ['frame'],
}

def _collect_particles(molecules, molecule_type, target_range):
    positions = []
    molecule_ids = []
    counts = []
    molecule_id = 0
    for molecule in molecules:
        if molecule_type != "" and molecule['molecule_type'] != molecule_type:
            continue
        n_in_mol = 0
        for i, position in enumerate(molecule['position']):
            if target_range != [] and i not in target_range:
                continue
            positions.append(position)
            molecule_ids.append(molecule_id)
            n_in_mol += 1
        if n_in_mol > 0:
            counts.append(n_in_mol)
            molecule_id += 1
    return positions, molecule_ids, counts

def _pair_count(n1, counts_1, counts_2, same_set, mode):
    if not same_set:
        return n1 * sum(counts_2)
    if mode == 'inter':
        total = sum(counts_1)
        return sum(n_in_mol * (total - n_in_mol) for n_in_mol in counts_1)
    if mode == 'intra':
        return sum(n_in_mol * (n_in_mol - 1) for n_in_mol in counts_1)
    return n1 * (n1 - 1)

def func(frame, type1='', type2='', dr=0.1, r_max=None, mode='', target_molecule=None, **kwargs):
    prop = frame.properties(ignoring_image=False)
    molecules = prop['molecules']
    box_size = prop['size']
    if r_max is None:
        r_max = min(box_size) / 2

    if type1 == '' and target_molecule not in (None, ''):
        type1 = target_molecule
    if type1 == '' and molecules:
        type1 = molecules[0]['molecule_type']
    if type2 == '':
        type2 = type1

    type1_name, range1 = protein_range_split(type1)
    type2_name, range2 = protein_range_split(type2)
    positions_1, mol_ids_1, counts_1 = _collect_particles(molecules, type1_name, range1)
    same_set = type1_name == type2_name and range1 == range2
    if same_set:
        positions_2, mol_ids_2, counts_2 = positions_1, mol_ids_1, counts_1
    else:
        positions_2, mol_ids_2, counts_2 = _collect_particles(molecules, type2_name, range2)

    n_bins = max(int(math.ceil(r_max / dr)), 1)
    radius_list = (np.arange(n_bins) + 0.5) * dr
    empty = Vector(
        data=np.zeros(n_bins),
        title='RDF',
        x_label='r (nm)',
        y_label='g(r)',
        x_axis=radius_list
    )
    n1 = len(positions_1)
    n2 = len(positions_2)
    if n1 == 0 or n2 == 0:
        return empty

    exclude_mode = 0
    if same_set:
        if mode == '':
            exclude_mode = 1
        elif mode == 'inter':
            exclude_mode = 2
        elif mode == 'intra':
            exclude_mode = 3

    hist = np.zeros(n_bins, dtype=np.int32)
    hist_device = cuda.to_device(hist)

    @cuda.jit
    def count_pairs(pos1, pos2, mol1, mol2, box, n_bins_, dr_, exclude, hist_):
        i = cuda.grid(1)
        while i < pos1.shape[0]:
            p1 = pos1[i]
            for j in range(pos2.shape[0]):
                if exclude == 1 and i == j:
                    continue
                if exclude == 2 and mol1[i] == mol2[j]:
                    continue
                if exclude == 3 and mol1[i] != mol2[j]:
                    continue
                if exclude == 3 and i == j:
                    continue
                dx = abs(p1[0] - pos2[j][0])
                dy = abs(p1[1] - pos2[j][1])
                dz = abs(p1[2] - pos2[j][2])
                if dx > box[0] / 2:
                    dx = box[0] - dx
                if dy > box[1] / 2:
                    dy = box[1] - dy
                if dz > box[2] / 2:
                    dz = box[2] - dz
                distance = math.sqrt(dx * dx + dy * dy + dz * dz)
                if distance <= 0:
                    continue
                bin_index = int(distance / dr_)
                if bin_index < n_bins_:
                    cuda.atomic.add(hist_, bin_index, 1)
            i += cuda.gridsize(1)

    count_pairs[1024, 128](
        cuda.to_device(np.array(positions_1, dtype=np.float32)),
        cuda.to_device(np.array(positions_2, dtype=np.float32)),
        cuda.to_device(np.array(mol_ids_1, dtype=np.int32)),
        cuda.to_device(np.array(mol_ids_2, dtype=np.int32)),
        cuda.to_device(np.array(box_size, dtype=np.float32)),
        n_bins,
        float(dr),
        exclude_mode,
        hist_device
    )
    hist = hist_device.copy_to_host().astype(float)

    volume = box_size[0] * box_size[1] * box_size[2]
    n_pair = _pair_count(n1, counts_1, counts_2, same_set, mode)
    shell_volume = 4.0 / 3.0 * math.pi * (
        (radius_list + dr / 2) ** 3 - (radius_list - dr / 2) ** 3
    )
    density = n_pair / volume if volume > 0 and n_pair > 0 else 1.0
    rdf = hist / (density * shell_volume)

    return Vector(
        data=rdf,
        title='RDF',
        x_label='r (nm)',
        y_label='g(r)',
        x_axis=radius_list
    )
