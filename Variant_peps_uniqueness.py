import read_fasta_file
from itertools import islice
import protein_digestor
import os

infile = "M_avium_Mavium_hominissuis_variant_proteome_search_082021_PSMs.txt"
fasta = "Mavium_hominissuis_GCF_000007865.1_ASM786v1_protein.fasta"
v_fasta = "../Mavium_hominissuis/Mavium_hominissuis_variant_proteins_DB_071421.fasta"

inpath = '.'

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

for fasta in os.listdir(os.path.join(inpath)):
    if os.path.isfile(os.path.join(inpath, fasta)):
        if fasta.split('.')[-1] == 'fasta':
            gen_peps(os.path.join(inpath, fasta))            

output = []
get_header(infile)
with open(infile) as file:
    for i in islice(file, 1, None):
        split_i = i.rstrip().split('\t')
        peps = split_i[4].strip('"').split('.')[1].upper()
        if peps in digested_pep:
            if len(digested_pep[peps]) == 1:
                for pro in digested_pep[peps]:
                    output.append(split_i + [pro])

##for o in output:
##   print  ('\t'.join(o) + '\n')

outfile = "{}_Final.txt".format(infile.rstrip('.txt'))
with open(outfile, 'w') as outf:
    outf.writelines('\t'.join(i) + '\n' for i in header)
    outf.writelines('\t'.join(i) + '\n' for i in output)
