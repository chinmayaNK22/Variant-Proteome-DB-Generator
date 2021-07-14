import read_fasta_file
from itertools import islice
import argparse
from datetime import date, datetime

parser = argparse.ArgumentParser(description='Custom generate proteome databases for protein variants such as SNPs from a list of SNPs corresponding to a protein.', add_help=True)

parser.add_argument('in_fasta', metavar='-FA', type=str, nargs='+', help='Proteome database of interest')

parser.add_argument('protein_variants_list', metavar='-F', type=str, nargs='+', help='A .txt (text) file with a list of protein variants')

#parser.add_argument('column_number_variants', metavar='-CV', type=str, nargs='+', help='Columns number in the text (.txt) file, where the variant/amino acid to be changed with site information is available for each protein (ex: R25P)')

#parser.add_argument('column_number_protein_accession_number', metavar='-CP', type=str, nargs='+', help='Columns number in the text (.txt) file, where the protein accession of variant/amino acid to be changed with site information is available')

args = parser.parse_args()


today = date.today()
dt = today.strftime("%m%d%y")

now = datetime.now()
current_time = now.strftime("%H%M%S")

def get_header_idx(infile):
    header = []
    with open(infile) as file:
        for i in islice(file,0,1):
            header = i.rstrip().split('\t')
            #header.append(i.rstrip())
    return header

def extract_SNPs(snps_file, SNPs_column_number, Proteins_column_number):
    snp = {}
    with open(snps_file) as infile:
        for i in islice(infile, 1, None):
            split_i = i.rstrip().split('\t')
            SNP = split_i[SNPs_column_number]
            Protein = split_i[Proteins_column_number]
            #print (SNP, Protein)
            if Protein not in snp:
                snp[Protein] = [SNP]
            else:
                snp[Protein].append(SNP)

    return snp

def variant_DB_generator(fasta, snps_file):#, SNPs_column, Proteins_column):
    head = get_header_idx(snps_file)
    print('Index\tColumn Name\n')
    for idx, columns in enumerate(head):
        print(str(idx)  + '\t' + columns)
    SNPs_column = input('Set the column index number where SNPs (ex:R35P) is present: ')
    Proteins_column = input('Set the column index number in which Protein Accession is present: ')
    snp_dicts = extract_SNPs(snps_file, int(SNPs_column), int(Proteins_column))
    outfile = '{0}_variant_DB_'.format(args.in_fasta[0].rstrip('.fasta')) + dt + '.fasta'
    with open(outfile, 'w') as outf:
        for rows in read_fasta_file.read_fasta(fasta):
           head = rows[0]
           seq = rows[1]
           accession = head.split(' ')[0]
           if accession in snp_dicts:
               final_SNPs = snp_dicts[accession]
               for snp in final_SNPs:
                   if not ':' in snp:
                       aa = snp[0]
                       new_aa = snp[-1]
                       site = ''.join(i for i in snp if i.isdigit())
                       #print (snp, site, len(seq)) #seq[int(site)-1])
                       try:
                           if seq[int(site)-1] == aa:
                               #print (accession, seq, final_SNP, seq[int(site)-1], site)
                               new_seq = seq[:int(site)-1] + seq[int(site)-1].replace(aa,new_aa) + seq[int(site):]
                               #print (head, seq, len(seq), seq[:int(site)-1] + '\n' + seq[int(site)-1].replace(aa,new_aa)+ '\' + seq[int(site):], final_SNP)
                               new_accession = accession + '_' + snp
                               new_header = '>' + new_accession + '|' + head
                               outf.write(new_header + '\n' + new_seq + '\n')
                           else:
                               print ('ERROR: Amino acid ' + aa + ' not found at position ' + str(site) + ' in protein ' + accession)
                       except:
                           pass
                           #print ('The variation ' + snp + ' cannot be incorporated in the protein ' + accession)
                   else:
                       pass
                       #print ('The variation ' + snp + ' cannot be incorporated in the protein ' + accession)
                  
variant_DB = variant_DB_generator(args.in_fasta[0], args.protein_variants_list[0])

## Usage
## python Variant_DB_generator.py Mycobacterium_tuberculosis_H37Rv_proteins_v3.fasta MTB_SNPs_010121.txt 1 2
