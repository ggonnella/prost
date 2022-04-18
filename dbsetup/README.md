# ProstDB automated database setup

Providing an automated system for setting up the RDBMS and starting the database
on different systems is very complex and out of the scope of this software.

However, a Snakefile ``prostdb.snake`` is provided, which may work in setting
up the database and starting the database server on some systems.

For this to work, the ``dbrootpass`` key of the configuration file
(``config.yaml``) must be set. Then the ProstDB snakefile
can be started using:
```
snakemake -j -s prostdb.snake connect
```

This should create the database, the database user and start the server.
If everything worked correctly a SQL prompt should appear
(you can exit it using ``quit;``).
