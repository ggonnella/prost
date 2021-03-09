# ProstDB

## Configuration

Change the values in `config.yaml` to adapt them to the local environment.

## Relative paths

Paths relative to the `datadir` directory are defined in the `common.snake`
file, which is imported by all other snakefiles, and are made available
as members of a `path` namespace (e.g. `path.dbdir`).

## Creating the database

Running `prostdb.snake` will initialize and create the DB
if it does not exist yet. It also starts or restarts the server if needed.

## NCBI taxonomy

Taxdump is the official dump of the NCBI taxonomy database,
(`ftp://ftp.ncbi.nih.gov/pub/taxonomy`). It is a tar.gz archive containing
multiple files. At Nov.2019 the compressed size was 49 Mb, uncompressed 344 Mb.
For more information (what is in each file) see
`ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt`.

Running `download_ncbi_taxonomy.snake` will update the dump files of
NCBI taxonomy, if necessary (i.e. if new files are available on the FTP server).
Nothing is downloaded otherwise.

Running `prostdb_ncbi_taxonomy.snake` will re-create the NCBI taxonomy DB
(dropping the tables and re-inserting all data) if new dumps are available
(i.e. if new files were output by `update_taxdump.snake`).  Nothing will happen
otherwise. Requirement: the "prostdb" database must exist and the connection
socket available.

### Extracting a subtree

The following extracts a subtree under a given node (in this case 123455)
using a SQL recursive query:
```
snakemake -j -s prostdb_ncbi_taxonomy.snake get_subtree -C root=123455
```

## NCBI complete prokaryotic genomes

Automatically checks for updates (i.e. for new complete genomes)
from NCBI and downloads them.

```
snakemake -j -s update_genome_sequences.snake
```

### Accession tables

```
snakemake -j -s compute_accession_tables.snake
```

### Assembly summaries database

```
snakemake -j -s prost_db_assembly_summaries.snake
```

## Bacdive

This downloads the data. It is not an automatic system to check for
updates, since the download is very slow.

```
snakemake -j -s download_bacdive.snake -f all
```

# Update strategy

To update all data (except Bacdive) run the following:
```
snakemake -j -s download_genome_sequences.snake
snakemake -j -s compute_accession_tables.snake
snakemake -j -s prostdb.snake
snakemake -j -s prostdb_assembly_summaries.snake
snakemake -j -s download_ncbi_taxonomy.snake
snakemake -j -s prostdb_ncbi_taxonomy.snake
```
