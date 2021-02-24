#!/usr/bin/env python3
"""
Create all NCBI taxonomy database tables.

Usage:
  ./ncbi_taxonomy_db_init.py [options] <dbuser> <dbpass> <dbname> <dbsocket>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file

Options:
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
from sqlalchemy import create_engine
from ncbi_taxonomy_db import Base
from docopt import docopt
from schema import Schema, And
import os

def main(arguments):
  connstr = "".join(["mysql+mysqldb://", arguments["<dbuser>"],
                     ":", arguments["<dbpass>"], "@localhost/",
                     arguments["<dbname>"], "?unix_socket=",
                     arguments["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  Base.metadata.create_all(engine.engine)

def validated(arguments):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   str: object})
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(validated(arguments))
