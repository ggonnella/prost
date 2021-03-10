#!/usr/bin/env python3
"""
List values of a given stat from all taxa under a given node.

Usage:
  query_stat_for_subtree.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <root> <table> <statname>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  root:      tax id of the root
  table:     table where the stat is stored
  statname:  value of the statname column

Options:
  --up, -u X       go up X levels into the taxonomic tree
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""

from sqlalchemy import create_engine
from dbschema.ncbi_taxonomy_db import NtNode
from dbschema.ncbi_assembly_summary import NcbiAssemblySummary
from dbschema.stats import tablename2class
from docopt import docopt
from schema import Schema, And, Use, Or, Optional
import os
from sqlalchemy.orm import sessionmaker, aliased

def traverse_up(session, node, levels):
  """
  Return the taxid of the root of the subtree
  by traversing up <levels> levels from node <node>
  """
  while levels:
    levels -= 1
    node = session.query(NtNode.parent_tax_id).\
        filter(NtNode.tax_id==node).one()[0]
  return node

def subtree_taxids(session, node):
  """
  Return the taxid of all nodes under the given <node>
  """
  subtree = session.query(NtNode)
  subtree = subtree.filter(NtNode.tax_id==node)
  subtree = subtree.cte(name="subtree", recursive=True)
  parent = aliased(subtree, name="parent")
  children = aliased(NtNode, name="children")
  subtree = subtree.union(
      session.query(children).\
          filter(children.parent_tax_id == parent.c.tax_id))
  ids = session.query(NtNode.tax_id).select_entity_from(subtree).all()
  return [i[0] for i in ids]

def accessions_for_taxids(session, taxids):
  """
  Return the Refseq accessions of complete genomes for a list of taxids
  """
  accessions = session.query(NcbiAssemblySummary.accession).\
      filter(NcbiAssemblySummary.taxid.in_(taxids)).\
      filter(NcbiAssemblySummary.seqdb == "refseq").\
      filter(NcbiAssemblySummary.assembly_level == "Complete Genome").all()
  return [a[0] for a in accessions]

def values_for_accessions(session, accessions, statname, table):
  """
  Return all available values of a stat, given a list of accessions
  """
  klass = tablename2class[table]
  values = session.query(klass.value).\
      filter(klass.statname == statname).\
      filter(klass.accession.in_(accessions)).all()
  return [v[0] for v in values]

def get_session(args):
  connstr = "".join(["mysql+mysqldb://", args["<dbuser>"],
                     ":", args["<dbpass>"], "@localhost/",
                     args["<dbname>"], "?unix_socket=",
                     args["<dbsocket>"]])
  engine = create_engine(connstr)
  return sessionmaker(bind=engine)()

def main(args):
  session = get_session(args)
  root = traverse_up(session, args["<root>"], args["--up"])
  taxids = subtree_taxids(session, root)
  accessions = accessions_for_taxids(session, taxids)
  values = values_for_accessions(session, accessions, args["<statname>"],
                                 args["<table>"])
  for value in values:
    print(value)

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<root>": Use(int),
                   "<table>": And(str, len),
                   "<statname>": And(str, len),
                   "--up": Or(None, Use(int)),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<dbuser>": snakemake.config["dbuser"],
    "<dbpass>": snakemake.config["dbpass"],
    "<dbname>": snakemake.config["dbname"],
    "<dbsocket>": snakemake.input.socket,
    "<root>": snakemake.params.root,
    "<table>": snakemake.params.table,
    "<statname>": snakemake.params.statname,
    "--up": snakemake.params.get("up", None)}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
