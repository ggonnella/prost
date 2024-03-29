#!/usr/bin/env python3
"""
Compute the ID of the nodes at taxonomic ranks from a given node up

Usage:
  ncbi_taxonomy_ranks_for_id.py [options] {db_args_usage} <node>

Arguments:
{db_args}
  node:      tax id of the node for which to compute the ranks

Output:
  - IDs of kingdom, phylum, class, order, family, genus, species
  - on one line, tab separated
  - if <node> is over species level, for the missing nodes, "None" is output

Options:
{common}
"""

from schema import Use
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import snacli
from dbschema.ncbi_taxonomy_db import NtNode
from lib import db, scripts

Ranks = ["superkingdom", "phylum", "class", "order",
         "family", "genus", "species"]

def compute_ranks(args):
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"])
  Session = sessionmaker(bind=engine)
  session = Session()
  prev_node_id = None
  node_id = args["<node>"]
  rank_ids = {}
  while prev_node_id != node_id:
    node_q = session.query(NtNode).filter(NtNode.tax_id==node_id)
    node = node_q.one()
    node_id = node.tax_id
    #print(f"{node.tax_id} {node.rank}")
    rank_ids[node.rank] = node.tax_id
    prev_node_id = node_id
    node_id = node.parent_tax_id
  output = []
  for rank in Ranks:
    if rank in rank_ids:
      output.append(str(rank_ids[rank]))
    else:
      output.append("None")
  return output

def main(args):
  output = compute_ranks(args)
  print("\t".join(output))

def validated(args):
  return scripts.validate(args, db.args_schema, {"<node>": Use(int)})

with snacli.args(db.snake_args,
                 params = ["<node>"],
                 docvars={"common": scripts.args_doc,
                 "db_args": db.args_doc, "db_args_usage": db.args_usage},
                 version="0.1") as args:
  if args: main(validated(args))
