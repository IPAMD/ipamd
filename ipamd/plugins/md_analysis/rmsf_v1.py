"""
plugin for calculating RMSF of each bead along a trajectory
"""
import numpy as np
from ipamd.public.models.md import Molecule
from ipamd.public.models.data import Vector
from ipamd.public.utils.output import warning
from ipamd.public.utils.parser import frame_index_list, protein_range_split

def _align_positions(positions, ref_positions):
    rotation = Molecule.align(
        {'position': positions},
        {'position': ref_positions}
    )[1]
    current = np.asarray(positions, dtype=float)
    reference = np.asarray(ref_positions, dtype=float)
    center_current = np.mean(current, axis=0)
    center_reference = np.mean(reference, axis=0)
    return (current - center_current) @ rotation + center_reference

def _selected_positions(prop, n=None, target_molecule=None):
    molecules = prop['molecules']
    if n is not None:
        return [np.array(molecules[n]['position'], dtype=float)]

    target_molecule_type, target_range = protein_range_split(target_molecule)
    selected = []
    for molecule in molecules:
        if target_molecule_type != "" and molecule['molecule_type'] != target_molecule_type:
            continue
        positions = molecule['position']
        if target_range:
            positions = [positions[i] for i in target_range]
        selected.append(np.array(positions, dtype=float))
    return selected

def func(box, target_frame, n=None, target_molecule=None, **kwargs):
    indices = frame_index_list(target_frame)
    if len(indices) < 2:
        warning('RMSF needs at least two frames')
        return Vector(
            data=[],
            title='RMSF',
            x_label='Residue Index',
            y_label='RMSF (nm)'
        )

    frames = [box.frame(i).current_frame() for i in indices]
    n_molecules = len(frames[0].molecules)
    if n is None and (target_molecule is None or target_molecule == '') and n_molecules > 1:
        n = 0
        warning('More than one molecule in the frame, only the first molecule will be used for RMSF calculation')

    ref_prop = frames[0].properties(ignoring_image=True, filter=None if n is None else str(n))
    ref_list = _selected_positions(ref_prop, n=0 if n is not None else None, target_molecule=target_molecule)
    n_beads = ref_list[0].shape[0] if ref_list else 0
    total_error = np.zeros(n_beads, dtype=float)
    n_samples = 0

    for frame in frames[1:]:
        prop = frame.properties(ignoring_image=True, filter=None if n is None else str(n))
        current_list = _selected_positions(prop, n=0 if n is not None else None, target_molecule=target_molecule)
        for current, reference in zip(current_list, ref_list):
            aligned = _align_positions(current, reference)
            delta = aligned - reference
            total_error += np.sum(delta ** 2, axis=1)
            n_samples += 1

    rmsf = np.sqrt(total_error / n_samples) if n_samples > 0 else np.zeros(n_beads)
    return Vector(
        data=rmsf,
        title='RMSF',
        x_label='Residue Index',
        y_label='RMSF (nm)',
        x_axis=list(range(n_beads))
    )
