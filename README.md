# Variant-Proteome-DB-Generator
A python command line based script for generating proteome databases for variants proteins from a list of variants.

## Usage
>python Variant_DB_generator.py in_fasta protein_variants_list column_number_variants column_number_protein_accession_number

1. in_fasta - Proteome database of interest in .fasta format
2. protein_variants_list - A .txt (text) file with a list of protein variants
3. column_number_variants - Columns number in the text (.txt) file, where the variant/amino acid to be changed with site information is available for each protein (ex: R25P)
4. column_number_protein_accession_number - Columns number in the text (.txt) file, where the protein accession of variant/amino acid to be changed with site information is available

## Example
>python Variant_DB_generator.py Mycobacterium_tuberculosis_H37Rv_proteins_v3.fasta MTB_SNPs_010121.txt 1 2
