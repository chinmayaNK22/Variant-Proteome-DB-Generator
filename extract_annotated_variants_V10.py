import os
from itertools import islice
import argparse

parser = argparse.ArgumentParser(description='''Extract protein coding variants from snpEff annotated .vcf file and save it in .txt format''')

parser.add_argument('infiles', metavar='-i', type=str, nargs='+', help='snpEff annotated .vcf files')

parser.add_argument('-f', '--feature_table', type=str, nargs='+', help='If snpEff provides only locus tag info, we need to extract protein accession details from feature table downloaded from RefSeq ftp path')

args = parser.parse_args()

aa = {'*':'*','GLY': 'G', 'ALA': 'A', 'SER': 'S', 'PRO': 'P', 'VAL': 'V', 'THR': 'T', 'CYS': 'C', 'LEU': 'L',
      'ILE': 'I', 'ASN': 'N', 'ASP': 'D', 'GLN': 'Q', 'LYS': 'K', 'GLU': 'E', 'MET': 'M', 'HIS': 'H','PHE': 'F',
      'SEC': 'U', 'ARG': 'R', 'TYR': 'Y', 'TRP': 'W', 'PYL': 'O', 'TER':'Ter'}

feature = {}
def extract_feature_table(f_table):
    with open(f_table) as file:
        for i in islice(file,1,None):
            if i.rstrip().startswith('CDS'):
                split_i = i.rstrip().split('\t')
                locus_tag = split_i[-3]
                pro_accession = split_i[-9]
                if locus_tag not in feature:
                    feature[locus_tag] = [pro_accession]
                else:
                    feature[locus_tag].append(pro_accession)

output = {}
variant_type = {}
def extract_ann_info(vcf_list, raw_file):
    variant_locus = {}
    if "#CHROM" not in vcf_list:
        chrom = vcf_list[0]
        chrom_pos = vcf_list[1]
        ref_nucleo = vcf_list[3]
        alt_nucleos = vcf_list[4]
        alt_nucleo = alt_nucleos.split(',')
        ann_info = vcf_list[7]
        annotations = []
        if ann_info.split(';')[-1].startswith('LOF='):
            annotations  = ann_info.split(';')[-2].split(',')
        else:
            annotations  = ann_info.split(';')[-1].split(',')

        for alt in range(len(alt_nucleos.split(','))):
            a = annotations[alt]
            if a.split('|')[-6].startswith('p.') and len(a.split('|')[-6]) != 0:
                try:
                    v_type = a.split('|')[1]
                    v_impact = a.split('|')[2]
                    if not 'synonymous_variant' in v_type:
                        if v_type not in variant_type:
                            variant_type[v_type] = [a.split('|')[3] + '@' + a.split('|')[-6]]
                        else:
                            variant_type[v_type].append(a.split('|')[3] + '@' + a.split('|')[-6])

                        if a.split('|')[3] not in variant_locus:
                            variant_locus[a.split('|')[3]] = [raw_file + '@' + chrom + '@' + chrom_pos + '@' + ref_nucleo + '@' + alt_nucleo[alt] + '@' + v_impact + '@' + v_type + '@' + a.split('|')[-6]]
                        else:
                            variant_locus[a.split('|')[3]].append(raw_file + '@' + chrom + '@' + chrom_pos + '@' + ref_nucleo + '@' + alt_nucleo[alt] + '@' + v_impact + '@' + v_type + '@' + a.split('|')[-6]) 
                    
                except:
                    print ('ERROR: ' + a + '    ' +  str(type(a)) + '    ' +  str(len(annotations)))

    for k, v in variant_locus.items():
        if k in feature:
            #protein = ''.join(pro for pro in feature[k])
            for j in v:
                try:
                    output_variants = ''
                    variants = j.split('@')[-1].lstrip('p.')
                    if 'ins' in variants:
                        ref_aa = ''
                        pos = ''
                        alt_aa_f = ''
                        for snp in range(len(variants.split('_'))):
                            ref_aa = ref_aa + aa[variants.split('_')[snp][0:3].upper()]
                            pos = pos + ':' + ''.join(word for word in variants.split('_')[snp] if word.isdigit())
                            alt_aa = variants.split('_')[snp].split(pos)[-1].split('ins')[-1]
                            alt_aa_f = ''.join(aa[alt_aa[i:i+3].upper()] for i in range(0,len(alt_aa),3) if alt_aa[i:i+3].upper() in aa)
                            #print (variants, variants.split('_')[snp]," ",ref_aa1, pos, alt_aa)
                        output_variants = ref_aa.lstrip() + '-' + pos.lstrip(':') + '-' + alt_aa_f
                        
                    elif 'del' in variants:
                        ref_aa = ''
                        pos = ''
                        alt_aa_f = ''
                        for snp in range(len(variants.split('_'))):
                            ref_aa = ref_aa + aa[variants.split('_')[snp][0:3].upper()]
                            pos = pos + ':' + ''.join(word for word in variants.split('_')[snp] if word.isdigit())
                            alt_aa = variants.split('_')[snp].split(pos)[-1].split('del')[-1].upper()
                            alt_aa_f = ''.join(aa[alt_aa[i:i+3].upper()] for i in range(0,len(alt_aa),3) if alt_aa[i:i+3].upper() in aa)
                            #print (variants, variants.split('_')[snp]," ",ref_aa1, pos, alt_aa)
                        output_variants = ref_aa.lstrip() + '-' + pos.lstrip(':') + '-' + alt_aa_f
                        
                    else:
                        ref_aa = ''
                        pos = ''
                        alt_aa = ''
                        if variants[0:3].upper() in aa:
                            ref_aa = aa[variants[0:3].upper()]
                        pos = ''.join(word for word in variants if word.isdigit())
                        if variants.split(pos)[-1].upper() in aa:
                            alt_aa = aa[variants.split(pos)[-1].upper()]

                        output_variants = ref_aa + '-' + pos + '-' + alt_aa
                    
                    op = j.split('@')[1] + '@' + j.split('@')[2] + '@' + j.split('@')[3] + '@' + j.split('@')[4] + '@' + j.split('@')[5] + '@' + k + '@' + ''.join(pro for pro in feature[k]) + '@' + j.split('@')[-2] + '@' + j.split('@')[-1] + '@' + output_variants
                    if op not in output:
                        output[op] = [j.split('@')[0]]
                    else:
                        output[op].append(j.split('@')[0])
        
                except:
                    print ("ERROR: Couldn't process " + variants + " from " + j.split('@')[0])

def extract_variants(vcf_file):
    c = 0
    if vcf_file.split('.')[-1] == 'vcf':
        with open(vcf_file) as file:
            for i in file:
                split_i = i.rstrip().split('\t')
                if len(split_i) == 10:
                    if split_i[4] != '.':
                        raw_read = os.path.split(vcf_file)[-1].rstrip('_ann.vcf')
                        extract_ann_info(split_i, raw_read)
                        #output.append(split_i)
                        c += 1
        print ('From file ' + vcf_file + ' ' + str(c) + ' variants were extracted')


    elif vcf_file.split('.')[-2] + '.' + vcf_file.split('.')[-1] == 'vcf.gz':
        for i in gzip.open(vcf_file, "rt"):
            split_i = i.rstrip().split('\t')
            if len(split_i) == 10:
                if split_i[4] != '.':
                    raw_read = os.path.split(vcf_file)[-1].rstrip('_ann.vcf.gz')
                    extract_ann_info(split_i, raw_read)
                    #output.append(split_i)
                    c += 1
        print ('From file ' + vcf_file + ' ' + str(c) + ' variants were extracted')

def extract_ann_vcf(inpath, feature_table):
    extract_feature_table(feature_table)
    if os.path.isfile(os.path.join(inpath)):
        if inpath.split('.')[-1] == 'vcf' or inpath.split('.')[-2] + '.' + inpath.split('.')[-1] == 'vcf.gz':
            extract_variants(os.path.join(inpath))
        else:
            print ('ERROR: The input file is not a vcf file')
    else:
        for files in os.listdir(os.path.join(inpath)):
            if os.path.isfile(os.path.join(inpath,files)):
                if files.split('.')[-1] == 'vcf' or files.split('.')[-2] + '.' + files.split('.')[-1] == 'vcf.gz':
                    extract_variants(os.path.join(inpath, files))
                else:
                    print ('ERROR: The input file is not a vcf file')

extract_ann_vcf(args.infiles[0], args.feature_table[0])

print (len(variant_type))

for k, v in variant_type.items():
        print (k, len(v))

results = []
for k, v in output.items():
    raw_files = ';'.join(r for r in v)
    results.append([raw_files] + [str(len(v))] + k.split('@'))
          
outfile = 'SnpEff_annotated_variants.txt'
with open(outfile, 'w') as outf:
    outf.write('File\tNo. of Files\tChromosome\tChromosome Position\tReference Nucleotide\tAlternative Nucleotide\tAnnotation Impact\tLocus_tag\tProtein_Accession\tVariant_Type\tSNP(Three Letter)\tSNP(Single Letter)\n')
    outf.writelines('\t'.join(i) + '\n' for i in results)
