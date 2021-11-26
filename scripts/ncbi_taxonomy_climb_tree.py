#!/usr/bin/env python3
"""
Compute the ID from a given node up to the root

Usage:
  ncbi_taxonomy_climb_tree.py [options] {db_args_usage} <node>

Arguments:
{db_args}
  node:      tax id of the node for which to start

Options:
{common}
"""

from sqlalchemy import create_engine
from dbschema.ncbi_taxonomy_db import NtNode
from docopt import docopt
from schema import Use
from lib import snake, db, scripts
from sqlalchemy.orm import sessionmaker, aliased

def main(args):
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"])
  Session = sessionmaker(bind=engine)
  session = Session()
  prev_node_id = None
  node_id = args["<node>"]
  while prev_node_id != node_id:
    node_q = session.query(NtNode).filter(NtNode.tax_id==node_id)
    node = node_q.one()
    node_id = node.tax_id
    print(node_id)
    prev_node_id = node_id
    node_id = node.parent_tax_id

def validated(args):
  return scripts.validate(args, db.args_schema, {"<node>": Use(int)})

if "snakemake" in globals():
  args = snake.args(snakemake, db.snake_args, params = ["<node>"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(db_args = db.args_doc, common = scripts.args_doc,
    db_args_usage = db.args_usage), version="0.1")
  main(validated(args))
