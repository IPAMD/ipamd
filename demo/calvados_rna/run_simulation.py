from ipamd import App
from ipamd.public.models.sequence import RNASequence
from ipamd.public.models.md import Unit, Environment

app = App('calvados_rna', gpu_id=3).use('Calvados2')
app.use('Calvados_rna', override=False)
box = app.builder.new_simulation_box(30, 30, 30)
seq = RNASequence('polyU', 'U' * 30)
mol_rna = app.builder.calvados_rna_from_sequence(seq, lambda index: (0, 0, index * 0.59))

density = 20
env = Environment()
env.set_ionic_strength(density * 1e-3)
env.set_temperature(293)
box.in_solvent(env)
box.place_molecule_randomly(mol_rna, 1)
simulation = app.simulation.new_simulation(box, f'simulation_{density}_mmol_L', total_time=20*Unit.TimeScale.ns, snap_shot=200)
simulation.run()


rg = app.mdanalysis.batch_compute(
    'rg_v1',
    box=box,
    target_frame='51-200'
)
print(rg.data)