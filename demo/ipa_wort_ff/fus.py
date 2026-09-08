from ipamd import App
from ipamd.public.models.md import Unit
app = App('ipa_wort_ff', gpu_id=0).use('ipa_wort')
app.builder.load_example_molecule('fus')
box = app.builder.new_simulation_box(40, 40, 40).in_solvent('pure water')
mol = app.builder.camc_protein_from_pdb('fus.pdb')
box.place_molecule_periodically(mol, 1, 1, 1)
box.to_xml('init')
