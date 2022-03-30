This document illustrates the first step in order to use Prost
in a new system.

These consists in the database creation, download of the data
and loading of the data into the database.

# Database creation

The database is created and initalized using the following command:
```
snakemake -j -s prostdb.snake
```

# Data download

The following commands download the data from external sources:
```
# NCBI Taxonomy:
#     Size: (Values at Jan 22) download size: 55 Mb
#                              uncompressed size: 388 Mb (+ compressed)
#     Location: <prostdata>/ncbi_taxonomy
#     More info: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt
#
snakemake -j -s download_ncbi_taxonomy.snake

# NCBI Assemblies, complete genomes:
#
snakemake -j -s download_ncbi_assemblies.snake
```

# Loading the data into the database


```
# NCBI Taxonomy:
snakemake -j -s prostdb_ncbi_taxonomy.snake

# Compute accession tables:
# (this takes e.g. 8 hours, since it downloads and temporarily decompress
#  all genomic sequences, to extract the sequence accessions)
snakemake -j -s compute_accession_tables.snake

# Load the assembly summaries into the database
snakemake -j -s prostdb_assembly_summary.snake
```

# Bacdive

This downloads the data. It is not an automatic system to check for
updates, since the download is very slow.

For this to work, the username and password for Bacdive must be
added to the ``config.yaml`` file.

```
snakemake -j -s download_bacdive.snake -f all
```

