#!/usr/bin/env python3
"""
Create NCBI taxonomy database indices for a table

Usage:
  ./ncbi_taxonomy_create_indices.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <table>

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

from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Sequence, Column, Integer, String, Index, DateTime
from sqlalchemy.orm import relationship
import ncbi_taxonomy_db
from sqlalchemy import create_engine
from ncbi_taxonomy_db import NtName
from docopt import docopt
from schema import Schema, And, Use
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
                     n in ncbi_taxonomy_db.tablename2class,
                     Use(lambda n: ncbi_taxonomy_db.tablename2class[n])),
                   str: object})
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(validated(arguments))
