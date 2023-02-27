# Prost

Prost is a system for the computation, storage and comparison of genomic data.
The name derives from `PROkaryote genome STatistics`.

The core component of the system is a database, named ProstDB, which stores:

- metadata about genome assemblies
- attributes: any kind of information collected or computed,
              about all or some of the genome assemblies
- computation tracing information
- taxonomic data (a mirror of NCBI taxonomy)
- phenotype data

ProstDB is implemented using the MariaDB database management system.

Furthermore, the system consists in programs, which allow to create, update,
compute new attributes, interact with the database and analyse the data
contained in it.

These are organized as follows:

- ProstPy
  a library of Python scripts, with a command line interface, as well
  as a SnakeMake API, which support the operations on the ProstDB database and
  on genomic data

- ProstPlugins
  plugins, written in Python, Nim and Rust that compute or collect data
  about the assemblies (attributes); each plugin computes one or multiple
  attributes; the user can write own plugins, according to the provided
  documentation (``docs/prost_plugins.md``); scripts for automatically checking
  the plugin interface are provided

- ProstSnakes
  a collection of pipeline scripts, based on SnakeMake, which
  automatically download and keep uptodata data from external databases, setup
  ProstDB and store the downloaded data in it, perform attributes computations
  and store the data in ProstDB, perform comparison tasks based on the ProstPy
  library and more.

A number of features has been created in separate packages
(installable using pip):

- multiplug: a plugin system which allows to write plugins for Python programs
             in Nim, Python, Rust and Bash; this is the base of ProstPlugins
- snacli: a system for creating scripts which can be run both from the command
          line (interactively) and from inside a Snakemake file;
          this is the base of the ProstPy interactive scripts
- attrtables: a dynamic system for storing entity attributes spread across
              multiple tables, including computation tracing metadata;
              this is the base of the ProstDB attributes storage
- prenacs: system for batch computing of entity attributes
           (using plugins based on ``multiplug``) and results storage
           (based on ``attrtables``) with provenance tracking;
           the library is used for the computation of genome assembly
           attributes and to store the results in the ProstDB
- findconfig: look for a configuration file in several possible positions
              on the file system

## Code organization

- ProstPy executable scripts are included in ``scripts``
- code common to multiple ProstPy scripts is stored under ``scripts/lib``
- ProstDB schema modules and further code to interact with the database
  is stored under ``scripts/dbschema``
- ProstPlugins are stored under ``plugins``
- ProstSnakes are stored under ``snakes``

## Documentation

For further information, see the documentation under ``docs``:
- ``setup.md`` describes how to install and setup Prost on a new system,
  including how to create the database, download the data and load it into the
  database
- ``data_update.md`` describes how to automatically update the data from
  external sources
- ``prost_plugins.md`` describes how to implement a plugin for attribute
  computation

## Acknowledgements
This software has been originally created in context of the DFG project GO 3192/1-1
“Automated characterization of microbial genomes and metagenomes by collection and verification of association rules”.
The funders had no role in study design, data collection and analysis.
