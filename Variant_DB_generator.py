import read_fasta_file
from itertools import islice
import argparse
from datetime import date, datetime

parser = argparse.ArgumentParser(description='Custom generate proteome databases for protein variants such as SNPs from a list of SNPs corresponding to a protein.', add_help=True)

parser.add_argument('in_fasta', metavar='-FA', type=str, nargs='+', help='Proteome database of interest')

parser.add_argument('protein_variants_list', metavar='-F', type=str, nargs='+', help='A .txt (text) file with a list of protein variants')

parser.add_argument('column_number_variants', metavar='-CV', type=str, nargs='+', help='Columns number in the text (.txt) file, where the variant/amino acid to be changed with site information is available for each protein (ex: R25P)')

parser.add_argument('column_number_protein_accession_number', metavar='-CP', type=str, nargs='+', help='Columns number in the text (.txt) file, where the protein accession of variant/amino acid to be changed with site information is available')

args = parser.parse_args()


today = date.today()
dt = today.strftime("%m%d%y")

now = datetime.now()
current_time = now.strftime("%H%M%S")


snp = {}
def extract_SNPs(snps_file, SNPs_column_number, Proteins_column_number):    
    with open(snps_file) as infile:
        for i in islice(infile, 1, None):
            split_i = i.split('\t')
            SNP = split_i[SNPs_column_number]
            Protein = split_i[Proteins_column_number]
            #print (SNP, Protein)
            snp[Protein] = SNP

    return snp

def variant_DB_generator(fasta, snps_file, SNPs_column, Proteins_column):
    ##print (fasta, snps_file, SNPs_column, Proteins_column)
    snp_dicts = extract_SNPs(snps_file, int(SNPs_column)-1, int(Proteins_column)-1)
    outfile = '{0}_variant_DB_'.format(args.in_fasta[0].rstrip('.fasta')) + dt + '.fasta'
    with open(outfile, 'w') as outf:
        ##write1 = open('variant_DB.fasta', 'w')
        for rows in read_fasta_file.read_fasta(fasta):
           head = rows[0]
           seq = rows[1]
           accession = head.split('|')[1]
           if accession in snp_dicts:
               final_SNP = snp_dicts[accession]
               aa = final_SNP[0]
               new_aa = final_SNP[-1]
               site = ''.join(i for i in final_SNP if i.isdigit())
               #print (final_SNP, seq[int(site)-1])
               if seq[int(site)-1] == aa:
                   #print (accession, seq, final_SNP, seq[int(site)-1], site)
                   new_seq = seq[:int(site)-1] + seq[int(site)-1].replace(aa,new_aa) + seq[int(site):]
                   #print (head, seq, len(seq), seq[:int(site)-1] + '\n' + seq[int(site)-1].replace(aa,new_aa)+ '\' + seq[int(site):], final_SNP)
                   new_accession = accession + '_' + final_SNP
                   new_header = '>' + new_accession + '|' + head
                   outf.write(new_header + '\n' + new_seq + '\n')

##    write1.close()
                  
variant_DB = variant_DB_generator(args.in_fasta[0], args.protein_variants_list[0], args.column_number_variants[0], args.column_number_protein_accession_number[0])

##if __name__== "__main__":
##    variant_DB_generator(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

##fa = 'Mycobacterium_tuberculosis_H37Rv_proteins_v3.fasta'
##fl = 'MTB_SNPs_010121.txt'
##c1 = 1
##p1 = 2
##
##
##variant_DB_generator(fa,fl,c1,p1)

## Usage
## python Variant_DB_generator.py Mycobacterium_tuberculosis_H37Rv_proteins_v3.fasta MTB_SNPs_010121.txt 1 2
