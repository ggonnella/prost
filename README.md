# Prost

Prost is a system for the computation, storage and comparison of genomic data.
The name derives from `PROkaryote genome STatistics`.

The main components of the system are:

- ProstDB: a database, conceived using MariaDB (but it can be probably easily
  adapted to other SQL database management systems), for storing genomic
  metadata and the results of computations on genomes, including detailled
  computation tracing information

- ProstPy: a library of Python scripts, with a command line interface, as well
  as a SnakeMake API, which support the operations on the ProstDB database and
  on genomic data; the library supports plugins written in Python, Nim and Rust

- ProstSnakes: a collection of pipeline scripts, based on SnakeMake, which
  automatically download and keep uptodata data from external databases, store
  data in the ProstDB database, start computations and comparison tasks
  based on the ProstPy library etc

# Installation

## Configuration

The file `config.yaml` must be edited, and information must be entered,
such as the system user, which will operate the system and the directory where
to store the data.

Developer notes: If Prost is installed as a git repository: `config.yaml` is
contained in the repository, but without some of the necessary values
(`dataroot`, `sysuser`). The user shall set them up. To avoid that the changes
in this file are marked by git use the following:
```
git update-index --skip-worktree config.yaml
```

## Requirements

The following software is required:
- Python3
- several Python libraries (sqlalchemy, sqlalchemy-repr, icecream, mariadb,
gffutils, sh)
- MariaDB (it can also a non-root user installation)
- Nim and/or Rust for the Nim/Rust plugins

If MariaDB is installed without using the package manager, the path of the
shared library must be communicated to the Python connector; this can be
done as follows:
- if root: editing /etc/ld.so.conf and running ldconfig afterwards
- if non-root: adding the path of libmariadb.so to the `LD_LIBRARY_PATH`

## Data download

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

## ProstDB

The following command will initialize and create the ProstDB and load the data
into the database (except Bacdive, see below):
```
# Database initialization:
snakemake -j -s prostdb.snake

# NCBI Taxonomy:
snakemake -j -s prostdb_ncbi_taxonomy.snake

# Compute accession tables:
# (this takes e.g. 8 hours, since it downloads and temporarily decompress
#  all genomic sequences, to extract the sequence accessions)
snakemake -j -s compute_accession_tables.snake

# Load the assembly summaries into the database
snakemake -j -s prostdb_assembly_summary.snake
```

# Update strategy

To update all data (except Bacdive, see below) run the following:
```
# Re-download NCBI Taxonomy data if there is any update:
snakemake -j -s download_ncbi_taxonomy.snake
# Update the downloaded NCBI assemblies
snakemake -j -s download_ncbi_assemblies.snake
# Fire up the database server if needed:
snakemake -j -s prostdb.snake
# Load the data into the database (if it has been updated):
snakemake -j -s prostdb_ncbi_taxonomy.snake
# Update the accession tables
snakemake -j -s compute_accession_tables.snake
# Load the assembly summaries into the database
snakemake -j -s prostdb_assembly_summary.snake
```

# Developer notes

## Relative paths

Paths relative to the `datadir` directory are defined in the `common.snake`
file, which is imported by all other snakefiles, and are made available
as members of a `path` namespace (e.g. `path.dbdir`).

# Extracting a subtree

The following extracts a subtree under a given node (in this case 123455)
using a SQL recursive query:
```
snakemake -j -s prostdb_ncbi_taxonomy.snake get_subtree -C root=123455
```

# Bacdive

This downloads the data. It is not an automatic system to check for
updates, since the download is very slow.

```
snakemake -j -s download_bacdive.snake -f all
```

