import copy
import numpy as np
from ipamd import App
from ipamd.public.models.sequence import ProteinSequence

app = App('mdanalysis_operator_test', gpu_id=0).use('Calvados3')
seq = ProteinSequence('polyA', 'A' * 40)
mol = app.builder.linear_protein(seq)
box = app.builder.new_simulation_box(40, 40, 40).in_solvent('pure water')
box.place_molecule_at(mol, 0, 0, 0)

first_frame = box.current_frame()
for shift in range(1, 5):
    box.new_frame()
    for item in first_frame.molecules:
        box.current_frame().add_molecule(
            copy.deepcopy(item['prototype']),
            offset=[item['offset'][0] + shift, item['offset'][1], item['offset'][2]]
        )

momentum = app.mdanalysis.momentum_v1(box=box, target_frame=0)
assert list(momentum.meta['x_axis']) == ['px', 'py', 'pz']
assert np.allclose(momentum.data, [0, 0, 0])

temperature = app.mdanalysis.temperature_v1(box=box, target_frame=0)
assert temperature.data == 0

rdf = app.mdanalysis.rdf(box=box, target_frame=0, type1='polyA', dr=0.2, r_max=4, mode='intra')
assert rdf.data.shape[0] == 20
assert np.all(np.isfinite(rdf.data))

ripley = app.mdanalysis.ripley_v1(box=box, target_frame=0, start_d=1, end_d=3, l=True)
assert len(ripley.data) == 3

rmsd = app.mdanalysis.rmsd_v1(box=box, target_frame='0-4', n=0)
assert np.allclose(rmsd.data, 0, atol=1e-6)

rmsf = app.mdanalysis.rmsf_v1(box=box, target_frame='0-4', n=0)
assert np.allclose(rmsf.data, 0, atol=1e-6)

msd = app.mdanalysis.msd(box=box, target_frame='0-4', n=0)
assert np.allclose(msd.data, [lag ** 2 for lag in range(5)], atol=1e-5)

flory = app.mdanalysis.flory_monomer_v1(box=box, target_frame='0-4', target_molecule='polyA')
assert np.isfinite(flory.data)

print('operator migration tests passed')
print('temperature', float(temperature.data))
print('rdf mean', float(np.mean(rdf.data)))
print('msd', list(np.round(msd.data, 6)))
print('flory', float(flory.data))
