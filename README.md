# Variant-Proteome-DB-Generator
A python command line based script for generating proteome databases for variants proteins from a list of variants.

## How to use Variant-Proteome-DB-Generator
```
>python Variant_DB_generator.py in_fasta protein_variants_list column_number_variants column_number_protein_accession_number

usage: Variant_DB_generator.py [-h] -FA [-FA ...] -F [-F ...]

Custom generate proteome databases for protein variants such as SNPs from a
list of SNPs corresponding to a protein.

positional arguments:
  -FA         Proteome database of interest
  -F          A .txt (text) file with a list of protein variants

optional arguments:
  -h, --help  show this help message and exit
```
