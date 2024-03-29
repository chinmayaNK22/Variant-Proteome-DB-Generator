# Variant-Proteome-DB-Generator
A python command line based script for generating proteome databases for variants proteins from a list of variants.

## How to use Variant-Proteome-DB-Generator
```
>python Variant_DB_generator.py proteome_database.fasta protein_variants_list.txt

usage: Variant_DB_generator.py [-h] -FA [-FA ...] -F [-F ...]

Custom generate proteome databases for protein variants such as SNPs from a list of SNPs corresponding to a protein.

positional arguments:
  -FA         Proteome database of interest
  -F          A .txt (text) file with a list of protein variants

optional arguments:
  -h, --help  show this help message and exit
```

## How to extract the SnpEff annotated variants from vcf files
```
>python extract_annotated_variants.py test.vcf -f feature_table.txt

usage: extract_annotated_variants.py [-h]
                                     [-f FEATURE_TABLE [FEATURE_TABLE ...]]
                                     -i [-i ...]

Extract protein coding variants from snpEff annotated .vcf file and save it in .txt format

positional arguments:
  -i                    snpEff annotated .vcf files

optional arguments:
  -h, --help            show this help message and exit
  -f FEATURE_TABLE [FEATURE_TABLE ...], --feature_table FEATURE_TABLE [FEATURE_TABLE ...]
                        If snpEff provides only locus tag info, we need to
                        extract protein accession details from feature table
                        downloaded from RefSeq ftp path
```

## How to check the uniqueness of variant peptides from variant DB search
```
>python Variant_peps_uniqueness.py -h
usage: Variant_peps_uniqueness.py [-h]
                                  -ip [-ip ...] -rf [-rf ...] -vf [-vf ...]

Uniqueness of variant peptides identified from the DDA database search can be
checked by matching it with reference and variant proteome databases

positional arguments:
  -ip         Exported PSMs of variant database search from Proteome
              Discoverer
  -rf         Reference proteome database of the same species in fasta format
  -vf         Variant proteome database in fasta format used for the search

optional arguments:
  -h, --help  show this help message and exit
```
# Citation
Please cite these tools using the DOI [![DOI](https://zenodo.org/badge/326594796.svg)](https://zenodo.org/badge/latestdoi/326594796)
