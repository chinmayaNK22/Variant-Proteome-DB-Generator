import read_fasta_file
from itertools import islice
import protein_digestor
import os
import argparse

<<<<<<< HEAD
parser = argparse.ArgumentParser(description='''Uniqueness of variant peptides identified from the DDA database search can be checked by matching it with
                                reference and variant proteome databases''')

parser.add_argument('infile', metavar='-ip', type=str, nargs='+', help='Exported PSMs of variant database search from Proteome Discoverer')
parser.add_argument('ref_fasta', metavar='-rf', type=str, nargs='+', help='Reference proteome database of the same species in fasta format')
parser.add_argument('ver_fasta', metavar='-vf', type=str, nargs='+', help='Variant proteome database in fasta format used for the search')

args = parser.parse_args()


def get_header_index(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')
            if '"' in split_i[0]:
                try:
                    peps_idx = split_i.index('"Annotated Sequence"')
                    return peps_idx
                except:
                    peps_idx = split_i.index('"Sequence"')
                    return peps_idx
            else:
                try:
                    peps_idx = split_i.index("Annotated Sequence")
                    return peps_idx
                except:
                    peps_idx = split_i.index("Sequence")
                    return peps_idx
=======
infile = "Example/M_avium_Mavium_hominissuis_variant_proteome_search_082021_PSMs.txt"
fasta = "Example/Mavium_hominissuis_GCF_000007865.1_ASM786v1_protein.fasta"
v_fasta = "Example/Mavium_hominissuis/Mavium_hominissuis_variant_proteins_DB_071421.fasta"
>>>>>>> 13939f748b8d878cd6b79fd7b7aac3ebf35ffc1c


header = []
def get_header(infile):
    with open(infile) as file:
        for i in islice(file, 0,1):
            header.append(i.rstrip().split('\t') + ['Protein Accession'] + ['Peptide Type'])

def digestor(seq,iter_cleavage,min_length,max_length, enzyme):
    pep = ''
    if enzyme == 'Trypsin':
        pep = protein_digestor.Trypsin(seq,iter_cleavage,min_length,max_length)
    elif enzyme == 'Lysc':
        pep = protein_digestor.Lysc(seq,iter_cleavage,min_length,max_length)
    elif enzyme == 'Chymotrypsin':
        pep = protein_digestor.Chymotrypsin(seq,iter_cleavage,min_length,max_length)

    return pep

digested_pep = {}
def gen_peps(infasta):
    for rows in read_fasta_file.read_fasta(infasta):
        peps = [i for iter_cleavage in range(int(2) + 1) for i in digestor(rows[1],iter_cleavage,int(6),int(60),'Trypsin')]
        accession = ""
        if '|' in rows[0]:
            accession = rows[0].split('|')[0] + '\t' + 'V'
        else:
            accession = rows[0].split(' ')[0] + '\t' + 'N'
        for pep in peps:
            if pep not in digested_pep:
                digested_pep[pep] = [accession]

            else:
                digested_pep[pep].append(accession)

def check_variant_uniq(infile, ref_fasta, ver_fasta):
    gen_peps(os.path.join(ref_fasta))
    gen_peps(os.path.join(ver_fasta))

    a = get_header_index(infile)
    output = []
    get_header(infile)
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            peps = split_i[a].strip('"').split('.')[1].upper()
            if peps in digested_pep:
                if len(digested_pep[peps]) == 1:
                    for pro in digested_pep[peps]:
                        output.append(split_i + [pro])


    outfile = "{}_Final.txt".format(infile.rstrip('.txt'))
    with open(outfile, 'w') as outf:
        outf.writelines('\t'.join(i) + '\n' for i in header)
        outf.writelines('\t'.join(i) + '\n' for i in output)

if __name__== "__main__":
    check_variant_uniq(args.infile[0], args.ref_fasta[0], args.ver_fasta[0])
