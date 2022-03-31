This document describes how to automatically update the data contained in the
databases.

# Update strategy

To update all data run the following:
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

There is no automatic update for the data in Bacdive.

