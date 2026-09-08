"""
plugin for computing Ripley's K or L function
"""
import math
import numpy as np
from numba import cuda
from ipamd.public.models.data import Vector
from ipamd.public.utils.parser import protein_range_split

configure = {
    "schema": ['frame'],
}

def func(frame, start_d=1, step=1, end_d=10, l=False, ref=False, target_molecule=None, **kwargs):
    @cuda.jit
    def count_neighbors(pos_list, box_size, max_d, res_list):
        i = cuda.grid(1)
        n_atoms = pos_list.shape[0]
        while i < n_atoms:
            position_i = pos_list[i]
            neighbor_count = 0
            max_image_x = int(max_d * 2 / box_size[0]) + 1
            max_image_y = int(max_d * 2 / box_size[1]) + 1
            max_image_z = int(max_d * 2 / box_size[2]) + 1
            for j in range(n_atoms):
                if i == j:
                    continue
                position_j = pos_list[j]
                for image_x in range(-max_image_x, max_image_x + 1):
                    for image_y in range(-max_image_y, max_image_y + 1):
                        for image_z in range(-max_image_z, max_image_z + 1):
                            delta_x = position_i[0] - position_j[0] - image_x * box_size[0]
                            delta_y = position_i[1] - position_j[1] - image_y * box_size[1]
                            delta_z = position_i[2] - position_j[2] - image_z * box_size[2]
                            distance_square = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z
                            if distance_square < max_d * max_d:
                                neighbor_count += 1
            res_list[i] = neighbor_count
            i += cuda.gridsize(1)

    prop = frame.properties(ignoring_image=False)
    molecules = prop['molecules']
    box_size = prop['size']
    target_molecule_type, target_range = protein_range_split(target_molecule)

    position_list = []
    for molecule in molecules:
        if target_molecule_type != "" and molecule['molecule_type'] != target_molecule_type:
            continue
        for i, position in enumerate(molecule['position']):
            if target_range != [] and i not in target_range:
                continue
            position_list.append(position)

    n_atoms = len(position_list)
    distances = list(range(int(start_d), int(end_d) + 1, int(step)))
    values = []

    if n_atoms == 0:
        return Vector(
            data=[0.0] * len(distances),
            title="Ripley's L" if l else "Ripley's K",
            x_label='r (nm)',
            y_label="L(r) (nm)" if l else "K(r) (nm^3)",
            x_axis=distances
        )

    position_list_device = cuda.to_device(np.array(position_list, dtype=np.float32))
    box_size_device = cuda.to_device(np.array(box_size, dtype=np.float32))
    volume = box_size[0] * box_size[1] * box_size[2]

    for d in distances:
        if ref:
            if l:
                values.append(0.0)
            else:
                values.append(4.0 / 3.0 * d ** 3 * np.pi)
            continue
        res_list = np.zeros(n_atoms, dtype=np.float32)
        res_list_device = cuda.to_device(res_list)
        count_neighbors[1024, 128](position_list_device, box_size_device, float(d), res_list_device)
        res_list_device.copy_to_host(res_list)
        n = float(np.sum(res_list))
        k = n * volume / (n_atoms ** 2)
        if l:
            values.append(math.pow(k / math.pi / 4.0 * 3.0, 1.0 / 3.0) - d)
        else:
            values.append(k)

    return Vector(
        data=values,
        title="Ripley's L" if l else "Ripley's K",
        x_label='r (nm)',
        y_label="L(r) (nm)" if l else "K(r) (nm^3)",
        x_axis=distances
    )
