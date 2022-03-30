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

## Documentation

For further information, see the documentation under ``manuals``:
- ``setup.py`` describes how to install and setup Prost on a new system
- ``init.py`` describes how to create the database, download the data and
              load it into the database, after the setup is completed
- ``update.py`` describes how to automatically update the data from external
                sources
