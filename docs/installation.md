# Installing and setup Prost

The present document illustrates how to setup Prost for a new system.

In particular Python (and a number of Python libraries), as well as
a relational database management system (RDBMS) are necessary.

The system was developed and tested using MariaDB as RDBMS. However,
any other RDBMS system could probably be used, with some adjustments.

## Requirements

The following software is required:
- Python3
  - several Python libraries, which can be installed using
  ``pip install -r requirements.txt``
- MariaDB
  (it can also be a non-root user installation)

Optional:
- Nim for Nim plugins
  - Nimble library ``zip`` for the example plugin
    ``nimble install zip``
- Rust for Rust plugins

## MariaDB library setup

If MariaDB is installed without using the package manager, the path of the
shared library must be communicated to the Python connector.

This can be done in one of the two following ways:
- as root user:
  editing /etc/ld.so.conf and running ldconfig afterwards
- as non root user:
  adding the path of libmariadb.so to the `LD_LIBRARY_PATH`

## Database creation

After MariaDB installation:
- a database must be created (e.g. named ``prostdb``)
- a database user with privileges on that database must be created
- the server must be started

## Configuration

The file `config.yaml` in the main project directory must be edited,
and information must be entered,
such as the system user, which will operate the system and the directory where
to store the data.

The file contains help text which explains what to enter under each
configuration key.

### Developer note:

If Prost is installed as a git repository: `config.yaml` is
contained in the repository, but without some of the necessary values
(`dataroot`, `sysuser`). The user shall set them up. To avoid that the changes
in this file are marked by git use the following:
```
git update-index --skip-worktree config.yaml
```

### Test suite database

To use the test suite, a test database (e.g. ``prosttestdb``) is necessary too.
The same user for the main database is used, and it must be granted full
privileges on the test database as well.
