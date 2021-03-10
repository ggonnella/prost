#!/usr/bin/env python3
"""
Given a genome, compute a statistics and compare it to taxonomic groups

Usage:
  cmpstat_taxtree.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <table> <statname> <genome> <module> <taxid> <up>

Arguments:
  table: table where the statistics is stored
  statname: name of the statistics
  genome: fas.gz file of the genome of interest
  module: module to compute the statistics
  taxid: tax ID of the genome of interest
  up: number of tax levels to climb up

Options:
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""

from schema import Schema, And, Use, Optional
from docopt import docopt
from dbschema.ncbi_taxonomy_db import NtName
import query_stat_for_subtree as qsfs
import cmpstat_distribution as cmpd
import genomestat_for_fas_gz as gffg
import os

def get_taxname(session, taxid):
  return session.query(NtName.name_txt).\
      filter(NtName.tax_id==taxid).\
      filter(NtName.name_class=="scientific name").one()[0]

def compare_to_subtree(session, node, value, statname, table):
 taxids=qsfs.subtree_taxids(session, node)
 accessions=qsfs.accessions_for_taxids(session, taxids)
 values=qsfs.values_for_accessions(session, accessions, statname, table)
 values=[float(v) for v in values]
 taxname = get_taxname(session, node)
 print("=================================================")
 print(f"Comparison of {statname} to the group: {taxname} (taxid: {node})")
 print(f"Value for the genome of interest: {value}")
 cmpd.cmp_to_distri(value, values)

def main(args):
  session = qsfs.get_session(args)
  assert(session)
  node = args["<taxid>"]
  value = float(gffg.compute_value(args["<module>"], args["<genome>"]))
  for i in range(args["<up>"]):
    compare_to_subtree(session, node, value,
        args["<statname>"], args["<table>"])
    node = qsfs.traverse_up(session, node, 1)

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<taxid>": Use(int),
                   "<table>": And(str, len),
                   "<statname>": And(str, len),
                   "<module>": open,
                   "<genome>": open,
                   "<up>": Use(int),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<dbuser>": snakemake.config["dbuser"],
    "<dbpass>": snakemake.config["dbpass"],
    "<dbname>": snakemake.config["dbname"],
    "<dbsocket>": snakemake.input.socket,
    "<taxid>": snakemake.params.taxid,
    "<table>": snakemake.params.table,
    "<statname>": snakemake.params.statname,
    "<module>": snakemake.input.module,
    "<genome>": snakemake.input.genome,
    "<up>": snakemake.params.up}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
