import subprocess
import os
import re
from shutil import copyfile
from pybioseq.fasta import FastaWriter
from ipamd.public.models.sequence import ProteinSequence
from ipamd.public.utils.output import info

configure = {
    "resource": ['persistency_dir'],
}
def func(
    protein_sequence : ProteinSequence,
    persistency_dir : str
):
    fasta_path = os.path.join(persistency_dir, protein_sequence.name + '.fasta')
    output_dir = os.path.join(persistency_dir, protein_sequence.name + '-af')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    writer = FastaWriter(fasta_path)
    writer.write(protein_sequence.name,"",protein_sequence.sequence)
    info(f'running alphafold2 on sequence {protein_sequence.sequence}')
    subprocess.run(
        [
            'colabfold_batch',
            fasta_path,
            output_dir,
            '--amber'
        ],
        capture_output=True,
        check=True
    )
    info('Alphafold2 finished')
    files = os.listdir(output_dir)
    for file in files:
        pattern = re.compile(rf'^{protein_sequence.name}_relaxed.*001.*.pdb$')
        if pattern.match(file):
            pdb_path = os.path.join(output_dir, file)
            new_pdb_path = os.path.join(persistency_dir, f"{protein_sequence.name}.pdb")
            copyfile(pdb_path, new_pdb_path)
    return f"{protein_sequence.name}.pdb"
