from ipamd import App
from ipamd.public.models.md import Unit
app = App('example', gpu_id=6).use('Calvados3')
app.builder.load_example_molecule('fus')
box = app.builder.new_simulation_box(40, 40, 40).in_solvent('pure water')
mol = app.builder.protein_from_pdb('fus.pdb', rigid_range='285-371')
#if you want to use the elastic network, uncomment the following code:
#
#mol, extra_ff = app.builder.gen_elastic_network(mol, max_gap=4)
#app.use(extra_ff, override=False)
#
box.place_molecule_randomly(mol, 20)
simulation = app.simulation.new_simulation(
    box, 'simulation', total_time=10*Unit.TimeScale.ns, snap_shot=20, minimize_energy=False
)
simulation.run()
mda = app.mdanalysis
dp = app.data_process
dp.plot(
    mda.batch_compute(
        "contact_map_v1",
        target_frame='10-20',
        box=box,
        threshold=4,
        mode='inter',
    ),
    save_figure=True
)
dp.plot(
    mda.batch_compute(
        'rdf',
        box=box,
        target_frame='10-20',
        type1='fus',
        dr=0.2,
        r_max=20,
        mode='inter',
    ),
    save_figure=True
)
dp.plot(
    mda.msd(
        box=box,
        target_frame='1-20',
        target_molecule='fus',
        dt=0.5,
    ),
    save_figure=True
)
dp.plot(
    mda.batch_compute(
        'density_align',
        box=box,
        target_frame='10-20',
        target_molecule='fus',
        direction='X',
        d=0.5,
    ),
    save_figure=True
)
