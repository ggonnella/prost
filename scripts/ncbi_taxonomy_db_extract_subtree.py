#!/usr/bin/env python3
"""
Extract subtree of NCBI taxonomy under a node

Usage:
  ncbi_taxonomy_extract_subtree.py [options] {db_args_usage} <root>

Arguments:
{db_args}
  root:      tax id of the root

Options:
{common}
"""

from sqlalchemy import create_engine
from dbschema.ncbi_taxonomy_db import NtNode
from schema import Use
from lib import db, scripts
from sqlalchemy.orm import sessionmaker, aliased
import snacli

def main(args):
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"])
  Session = sessionmaker(bind=engine)
  session = Session()
  subtree = session.query(NtNode)
  subtree = subtree.filter(NtNode.tax_id==args["<root>"])
  subtree = subtree.cte(name="subtree", recursive=True)
  parent = aliased(subtree, name="parent")
  children = aliased(NtNode, name="children")
  subtree = subtree.union(
      session.query(children).\
          filter(children.parent_tax_id == parent.c.tax_id))
  result = session.query(NtNode).select_entity_from(subtree).all()
  for r in result:
    print(r)

def validated(args):
  return scripts.validate(args, db.args_schema, {"<root>": Use(int)})

with snacli.args(db.snake_args,
                 params = ["<root>"],
                 docvars={"common": scripts.args_doc,
                 "db_args": db.args_doc, "db_args_usage": db.args_usage},
                 version="0.1") as args:
  if args: main(validated(args))
