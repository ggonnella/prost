#!/usr/bin/env python3
"""
Given a genome, compute a statistics and compare it to taxonomic groups

Usage:
  seqattr_cmp_taxtree.py [options] {db_args_usage}
                         <attribute> <genome> <plugin> <taxid> <up>

Arguments:
{db_args}
  attribute:    name of the attribute
  genome:       fas.gz file of the genome of interest
  plugin:       plugin to compute the statistics
  taxid:        tax ID of the genome of interest
  up:           number of tax levels to climb up

Limitations:
  - (1) only for single-value attributes
  - (2) only for numeric (integer or float) attributes

Options:
{common}
"""

from schema import And, Use
from docopt import docopt
from dbschema.ncbi_taxonomy_db import NtName
import taxtree_query
import cmpstat_distribution
import batch_compute
from lib import snake, db, scripts, mod
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

def get_taxname(session, taxid):
  return session.query(NtName.name_txt).\
      filter(NtName.tax_id==taxid).\
      filter(NtName.name_class=="scientific name").one()[0]

def compare_to_subtree(session, engine, node, value, attribute):
 taxids=taxtree_query.subtree_taxids(session, node)
 accessions=taxtree_query.accessions_for_taxids(session, taxids)
 acc_values=taxtree_query.values_for_accessions(\
     session, engine, accessions, attribute)
 values=[float(v[0]) for v in acc_values.values()]
 taxname = get_taxname(session, node)
 print("=================================================")
 print(f"Comparison of {attribute} to the group: {taxname} (taxid: {node})")
 print(f"Value for the genome of interest: {value}")
 cmpstat_distribution.cmp_to_distri(value, values)

def main(args):
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"])
  session = sessionmaker(bind=engine)()
  plugin = mod.py_or_nim(args["<plugin>"], args["--verbose"])
  node = args["<taxid>"]
  results, logs = plugin.compute(args["<genome>"])
  value = float(results[plugin.OUTPUT.index(args["<attribute>"])])
  for i in range(args["<up>"]):
    compare_to_subtree(session, engine, node, value,
        args["<attribute>"])
    node = taxtree_query.traverse_up(session, node, 1)

def validated(args):
  return scripts.validate(args, db.args_schema, {"<taxid>": Use(int),
    "<attribute>": And(str, len), "<plugin>": open, "<genome>": open,
    "<up>": Use(int)})

if "snakemake" in globals():
  args = snake.args(snakemake, db.snake_args,
    params = ["<taxid>", "<attribute>", "<up>"],
    input = ["<plugin>", "<genome>"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(db_args = db.args_doc,
                db_args_usage = db.args_usage, common = scripts.args_doc),
                version="0.1")
  main(validated(args))
