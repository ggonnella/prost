## Update genome sequences

Automatically checks for updates (i.e. for new complete genomes)
from NCBI and downloads them.

```
snakemake --cores 4 -s update_genome_sequences.snake
snakemake --cores 4 -s update_accession_tables.snake
```

## Download Bacdive data

This downloads the data. It is not an automatic system to check for
updates, since the download is very slow.

```
snakemake --cores 4 -s download_bacdive.snake -f all
```
