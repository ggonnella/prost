#!/usr/bin/env python3
"""
Create all tables from a SqlAlchemy schema file.

Usage:
  db_create_tables.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <file>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  file:      python file containing the SqlAlchemy classes

Requirements:
  the file must define a global variable Base from declarative_base() e.g.
  Base = declarative_base(cls=PrettyRepresentableBase)

Options:
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
from sqlalchemy import create_engine
from docopt import docopt
from schema import Schema, And
import os
import importlib

def main(arguments):
  connstr = "".join(["mysql+mysqldb://", arguments["<dbuser>"],
                     ":", arguments["<dbpass>"], "@localhost/",
                     arguments["<dbname>"], "?unix_socket=",
                     arguments["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  spec = importlib.util.spec_from_file_location("models", arguments["<file>"])
  models = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(models)
  models.Base.metadata.create_all(engine.engine)

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
