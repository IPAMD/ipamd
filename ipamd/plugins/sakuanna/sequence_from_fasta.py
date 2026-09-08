import os.path
from pybioseq.fasta import FastaReader
from pybioseq.converter import standard_aa, standard_nmp, standard_dnmp
from ipamd.public.models.sequence import ProteinSequence, DNASequence, RNASequence
from ipamd.public.utils.output import error

configure = {
    "resource": ['persistency_dir']
}
def func(fasta_path, mol='protein', persistency_dir=None):
    fasta_full_path = os.path.join(persistency_dir, fasta_path)
    reader = FastaReader(fasta_full_path)
    if mol == 'protein':
        sequence_class = ProteinSequence
        check_function = standard_aa
    elif mol == 'dna':
        sequence_class = DNASequence
        check_function = standard_dnmp
    elif mol == 'rna':
        sequence_class = RNASequence
        check_function = standard_nmp
    else:
        raise ValueError(f"Unsupported molecule type: {mol}")

    sequence_list = []
    for name, _, sequence in reader.read():
        for char in sequence:
            if not check_function(char):
                error(f"Invalid character '{char}' found in sequence '{name}'.")
                raise ValueError()
        sequence_list.append(sequence_class(name, sequence))
    if len(sequence_list) == 1:
        return sequence_list[0]
    else:
        return tuple(sequence_list)
