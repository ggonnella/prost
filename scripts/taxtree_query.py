#!/usr/bin/env python3
"""
List values of a given attribute from all taxa under a given node.

Usage:
  taxtree_query.py [options] {db_args_usage}
                             <root> <attribute>

Arguments:
{db_args}
  root:       tax id of the subtree root
  attribute:  name of the attribute

Options:
  --up, -u X       go up X levels into the taxonomic tree
{common}
"""

from sqlalchemy import create_engine
from dbschema.ncbi_taxonomy_db import NtNode
from dbschema.ncbi_assembly_summary import NcbiAssemblySummary
from dbschema.attribute import AttributeValueTables
from lib import db, scripts
from schema import And, Use, Or
from sqlalchemy.orm import sessionmaker, aliased
import snacli

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

def values_for_accessions(session, engine, accessions, attribute):
  """
  Return all available values of an attribute, given a list of accessions
  """
  avt = AttributeValueTables(engine)
  t_sfx, vcolnames, ccolname, gcolname = avt.attribute_location(attribute)
  klass = avt.get_class(t_sfx)
  rows = session.query(klass).\
      filter(klass.c.accession.in_(accessions)).all()
  return {row.accession: [getattr(row, cn) for cn in vcolnames] for row in rows}

def main(args):
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"])
  session = sessionmaker(bind=engine)()
  root = traverse_up(session, args["<root>"], args["--up"])
  taxids = subtree_taxids(session, root)
  accessions = accessions_for_taxids(session, taxids)
  values = values_for_accessions(session, engine,
      accessions, args["<attribute>"])
  for acc, val in values.items():
    print(f"{acc}\t"+"\t".join([str(v) for v in val]))

def validated(args):
  return scripts.validate(args, db.args_schema, {"<root>": Use(int),
                   "<attribute>": And(str, len),
                   "--up": Or(None, Use(int))})

with snacli.args(db.snake_args,
           params = ["<root>", "<attribute>", "--up"],
                 docvars={"common": scripts.args_doc,
                 "db_args": db.args_doc, "db_args_usage": db.args_usage},
                 version="0.1") as args:
  main(validated(args))
