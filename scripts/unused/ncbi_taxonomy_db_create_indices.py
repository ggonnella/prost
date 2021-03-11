#!/usr/bin/env python3
"""
Create NCBI taxonomy database indices for a table

This was splitted from the table creation,
in order to defer it after the bulk import.

Usage:
  ncbi_taxonomy_create_indices.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <table>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  table:     table name

Options:
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""

from sqlalchemy import create_engine
from dbschema.ncbi_taxonomy_db import tablename2class
from docopt import docopt
from schema import Schema, And, Use, Optional
import os

def main(arguments):
  connstr = "".join(["mysql+mysqldb://", arguments["<dbuser>"],
                     ":", arguments["<dbpass>"], "@localhost/",
                     arguments["<dbname>"], "?unix_socket=",
                     arguments["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  arguments["<table>"].create_indices(engine)

def validated(arguments):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<table>": And(lambda n:
                     n in tablename2class,
                     Use(lambda n: tablename2class[n])),
                   Optional(str): object})
  return schema.validate(arguments)

if "snakemake" in globals():
  args = {
    "<dbuser>": snakemake.config["dbuser"],
    "<dbpass>": snakemake.config["dbpass"],
    "<dbname>": snakemake.config["dbname"],
    "<dbsocket>": snakemake.input.socket,
    "<table>": snakemake.params.table}
  main(validated(args))
elif __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(validated(arguments))
