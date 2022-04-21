# Prost setup guide

In order to use Prost for the first time, the following steps are necessary:

- Requirements setup:   installation of necessary software and environment setup
- Configuration:        a configuration file is prepared
- Database setup:       the database is initialized and prepared for use
- Data download:        external data is downloaded
- Metadata computation: some metadata is computed from the downloaded data
- Loading into the db:  the database is populated using downloaded data/metadata

## Requirements setup

Prost is written mainly in Python 3, which must be installed. Furthermore, some
Python libraries are necessary, which are listed in the requirements file
(``requirements.txt``) in the main directory of the package source code. They
can be e.g. installed using: ``pip install -r requirements.txt``.

Prost is based on plugins, which perform computations. The user can write
plugins in several programming languages: Python, Bash, Nim and Rust.
If Nim or Rust plugins are used, the necessary compilers must be installed.

A relational database management system (RDBMS) is necessary.
The system was developed and tested using MariaDB as RDBMS. However,
any other RDBMS system could probably be used, with some adjustments.
Note that the MariaDB can be installed also as non-root user.

### MariaDB library setup

If MariaDB is installed without using the package manager, the path of the
shared library must be communicated to the Python connector.

This can be done in one of the two following ways:
- as root user:
  editing /etc/ld.so.conf and running ldconfig afterwards
- as non root user:
  adding the path of libmariadb.so to the `LD_LIBRARY_PATH`

## Configuration

Prost requires a configuration file named ``prost.config.yaml``.
Prost looks for the file in the following paths, in the order
(and uses the first configuration file found):
- the path given in the ``XDG_CONFIG_DIR`` environment variable
- the directory ``.config`` under the user home directory
- the main directory of the Prost source code

The configuration consists of YAML file, representing a dictionary
with string keys and values.

An example configuration file is given in the source code directory
(``example.prost.config.yaml``). At minimum, the ``dataroot`` key
must be set (see below), then the file can be copied to one of the
locations indicated above, with the filename ``prost.config.yaml``.

If keys other than ``dataroot`` contain the string ``{dataroot}``,
this is replaced by the value of the ``dataroot`` entry.

### Data root

The ``dataroot`` key must be set to the path of the directory under which Prost
will store data downloaded from external sources, computed data, computation
state flags and (by default) the database data directory.

Developer note: the exact locations under ``dataroot`` of each kind of
information are determined in the ``common.snake`` module.

### Database connection configuration

It is possible to use an existing database server.
In this case the configuration file must contain the connection data:
- name of the database: `dbname`
- user and password of the database user: `dbuser`, `dbpass`
- hostname and port: `dbhost`, `dbport`
- socket file: `dbsocket`

If the database server does not exist yet, the desired connection data
must be entered. For example, one can use the values given in the example
configuration file. The instructions in the following section can then
be followed.

## Database setup

If a database server is already setup and ready for use with Prost, then
the connection data can be entered in the configuration file as explained
above and this section can be skipped.

For the database to be ready for use, the following steps are necessary:
- initialization of the data directory and creation of an administration
  account
- creation of a database and a database user to be used by Prost
- database server startup

These steps can be executed according the the instructions given in the
database management system manual (e.g. MariaDB) and the connection data can be
entered in the configuration file.

In alternative, the scripts and/or the snakefile ``prostdb.snake``, located in
the directory ``dbsetup`` can be used. These scripts are specific for MariaDB
and may not work on some systems.

### Database initialization

The database initialization creates the directory where the database data
is stored. The path to that directory must be added to the configuration
file under the key ``dbpath``. It shall not exist yet.

The initialization creates a database user with the same user name as the
operating system user name. The password to be used for this user
must be added to the configuration file under the key ``dbrootpass``.

The initialization is done using:
```
snakemake -j -s dbsetup/prostdb.snake initialize
```

### Database creation

The creation step follows the initialization and it creates the database to be
used for Prost and a user, which has full privileges on that database.

It requires the following variables in the configuration file (the
same as for the initialization step):
- the path to the database data directory (``dbpath``)
- the password of the database administrator user named as the operating
  system user (``dbrootpass``)

The step is executed using:
```
snakemake -j -s dbsetup/prostdb.snake create
```

### Database server startup and shutdown

The database server must be started, so that the connection is possible.
If the database has been created as explained above, the server can be started
using:
```
snakemake -j -s dbsetup/prostdb.snake start
```

A running database server can be stopped using:
```
snakemake -j -s dbsetup/prostdb.snake stop
```
The database server stop task requires the ``dbrootpass`` to be set in the
configuration.

### Test suite database

To run the test suite, a test database is necessary too. The database can
be created in the same data directory for the main database, using the
task ``create_test``. This requires the setup of a key ``testdbname``.
The same user (``dbuser``) as for the main database is used in the test
database, thus no further configuration entries are necessary.


This document illustrates the first step in order to use Prost
in a new system, i.e. the data download from external sources
and loading of the data into the database.

It assumes that a working relational database management system
has been installed and setup correctly (see ``installation.md``).

## Data download

The following commands download the data from external sources:
```
# NCBI Taxonomy:
#     Size: (Values at Jan 22) download size: 55 Mb
#                              uncompressed size: 388 Mb (+ compressed)
#     Location: <prostdata>/ncbi_taxonomy
#     More info: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt
#
snakemake -j -s snakes/download_ncbi_taxonomy.snake

# NCBI Assemblies, complete genomes:
#
snakemake -j -s snakes/download_ncbi_assemblies.snake
```

## Metadata computation

Sequence metadata which will be loaded into the database
must be computed for the genomic sequences, in particular the
sequence accessions corresponding to the assembly accessions.

For this task, each sequence file must be decompressed (temporarily).
The task may take several hours to complete (e.g. 8 hours).

```
snakemake -j -s snakes/compute_accession_tables.snake
```

### Loading the data into the database

Once the data has been downloaded and the metadata computed, the results
can be loaded into the database:

```
# NCBI Taxonomy:
snakemake -j -s snakes/prostdb_ncbi_taxonomy.snake

# Load the assembly summaries into the database
snakemake -j -s snakes/prostdb_assembly_summary.snake
```

# Bacdive (experimental)

Data from the Bacdive database can be downloaded. It is not an automatic system to check for
updates, since the download is very slow.

For this to work, the username (``bacdive_user``) and password
(``bacdive_pass``) for Bacdive must be added to the configuration file.

```
snakemake -j -s download_bacdive.snake -f all
```

