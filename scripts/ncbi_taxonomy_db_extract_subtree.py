#!/usr/bin/env python3
"""
Extract subtree of NCBI taxonomy under a node

Usage:
  ./ncbi_taxonomy_extract_subtree.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <root>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  root:      tax id of the root

Options:
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""

from sqlalchemy import create_engine
from ncbi_taxonomy_db import NtNode
from docopt import docopt
from schema import Schema, And, Use, Or
import os
from sqlalchemy.orm import sessionmaker, aliased

def main(arguments):
  connstr = "".join(["mysql+mysqldb://", arguments["<dbuser>"],
                     ":", arguments["<dbpass>"], "@localhost/",
                     arguments["<dbname>"], "?unix_socket=",
                     arguments["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  Session = sessionmaker(bind=engine)
  session = Session()
  subtree = session.query(NtNode)
  subtree = subtree.filter(NtNode.tax_id==arguments["<root>"])
  subtree = subtree.cte(name="subtree", recursive=True)
  parent = aliased(subtree, name="parent")
  children = aliased(NtNode, name="children")
  subtree = subtree.union(
      session.query(children).\
          filter(children.parent_tax_id == parent.c.tax_id))
  result = session.query(NtNode).select_entity_from(subtree).all()
  for r in result:
    print(r)

def validated(arguments):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<root>": Or(None, Use(int)),
                   str: object})
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(validated(arguments))
