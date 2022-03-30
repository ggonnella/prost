# Code organization

The tasks are defined in different snakemake files.

## ProstDB interaction

The ``prostdb.snake`` is used for the basic interaction with the database, i.e.
database creation, starting the server, stopping the server, connecting
interactively to the server.

## Download tasks

Files with the prefix ``download_..`` are used for tasks which download
information from external sources. Whenever possible, these shall only download
information which has changed, incrementally. This is not implemented for
Bacdive.

## Compute tasks

Files with prefix ``compute_..`` are used for computing data derived
from the data downloaded from external sources.

## Loading data into ProstDB

Files with the prefix ``prostdb_..`` are use to load the data from external
sources or from computations into the database.

## Attributes

The ``prostdb_attrs.snake`` is used for computing attributes and
loading the computation results into the database.

## Shared configuration

The ``config.yaml`` file contains configuration data. Some of these are preset,
some must be provided by the user. The ``common.snake`` file is used to load
and check the data contained in ``config.yaml`` and compute from it common
configuration variables.



