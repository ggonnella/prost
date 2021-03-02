## Update genome sequences

Automatically checks for updates (i.e. for new complete genomes)
from NCBI and downloads them.

```
snakemake -j -s update_genome_sequences.snake
snakemake -j -s update_accession_tables.snake
```

## Download Bacdive data

This downloads the data. It is not an automatic system to check for
updates, since the download is very slow.

```
snakemake -j -s download_bacdive.snake -f all
```

## Taxdump

Taxdump is the official dump of the NCBI taxonomy database,
(ftp://ftp.ncbi.nih.gov/pub/taxonomy). It is a tar.gz archive containing
multiple files. At Nov.2019 the compressed size was 49 Mb, uncompressed 344 Mb.
For more information (what is in each file) see
ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt

```
snakemake -j -s update_taxdump.snake
```

## Extract a subtree

```
snakemake -j -s update_ncbi_taxonomy_db.snake get_subtree -C root=123455
```
