#!/usr/bin/env python3
"""
Create all tables from a SqlAlchemy schema file.

Usage:
  db_create_tables.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <file> [<table>]

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  file:      python file containing the SqlAlchemy classes
  table:     create only specified table (optional, default: create all)

Requirements:
  the file must define a global variable Base from declarative_base() e.g.
  Base = declarative_base(cls=PrettyRepresentableBase)

Options:
  --drop, -d       drop table if exists (requires <table> argument)
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy import create_engine
from docopt import docopt
from schema import Schema, And, Optional, Or
import os
import importlib

def main(args):
  connstr = "".join(["mysql+mysqldb://", args["<dbuser>"],
                     ":", args["<dbpass>"], "@localhost/",
                     args["<dbname>"], "?unix_socket=",
                     args["<dbsocket>"]])
  engine = create_engine(connstr, echo=args["--verbose"])
  spec = importlib.util.spec_from_file_location("models", args["<file>"])
  models = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(models)
  if args["<table>"]:
    if args["--drop"]:
      query = "DROP TABLE IF EXISTS `{}`;".format(args["<table>"])
      Session = sessionmaker(bind=engine)
      session = Session()
      session.execute(query)
    models.Base.metadata.tables[args["<table>"]].create(bind=engine.engine)
  else:
    models.Base.metadata.create_all(engine.engine)

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<file>": open,
                   "<table>" : Or(None, str),
                   "--drop": Or(None, True, False),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
      "<dbuser>": snakemake.config["dbuser"],
      "<dbpass>": snakemake.config["dbpass"],
      "<dbname>": snakemake.config["dbname"],
      "<dbsocket>": snakemake.input.socket,
      "<file>": snakemake.input.schema,
      "<table>": snakemake.params.get("table", None),
      "--drop": snakemake.params.get("drop", None)
      }
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
