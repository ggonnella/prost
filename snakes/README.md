# ProstSnakes

ProstSnakes is a collection of Snakemake files (based on the ``snakemake`` tool),
which define automated tasks for interacting with the Prost
database, download external data, compute attributes and store the computation
results into the database. These operations are done by automating and
organizing calls to the ProstPy scripts (see ``scripts`` directory of the project).

The tasks are implemented in different modules (``*.snake`` files).

Besides automation of the tasks, the main design goal is incremental update.
So for example, only new assemblies are downloaded from NCBI Genomes each time.
Then, only attributes of the new downloaded genomes are computed.

## Configuration file and data location

The configuration file (``config.yaml`` in the main directory of the package)
defines the base directory where data
is stored (database, downloaded files, etc).  The subdirectories under this
point to be used for different purposes are defined in the common module
(``common.snake``), which is loaded by all other modules.

## Summary of the modules and tasks

Modules are called using ``snakemake -j -s <SNAKEFILE> [<TASK>]``.
If no ``<TASK>`` is specified, then the default task (``all``) is executed.

For brevity only the basename of the module (eg. without ``.snake`` extension)
and the task (if any) are mentioned in the table below:

| Module and task         | Purpose                                            |
|-------------------------|----------------------------------------------------|
| ``prostdb connect``     | interactively connect to the database      |
| ``download_ncbi_taxonomy`` | download the NCBI taxonomy database, if changed      |
| ``prostdb_ncbi_taxonomy``  | load downloaded NCBI taxonomy data into the database |
| ``download_ncbi_assemblies`` | incrementally download NCBI Genome assemblies      |
| ``compute_accession_tables`` | compute the accession tables from the downloaded assemblies |
| ``prostdb_assembly_summary`` | load NCBI assembly data to the database |
| ``prostdb_attrs``            | compute attributes for all assemblies |

