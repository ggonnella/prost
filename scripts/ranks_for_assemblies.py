#!/usr/bin/env python3
"""
Compute the taxonomy IDs for each rank for all accessions in a summary table.

Usage:
  rs_seq_accessions_for_assemblies.py [options] {db_args_usage} <summary>

Arguments:
{db_args}
  summary: NCBI assembly summary file

Options:
  --update, -u F            file with previous results, to update
{common}
"""

from sqlalchemy import create_engine, exc
from dbschema.ncbi_taxonomy_db import NtNode
from docopt import docopt
from schema import Or, And, Use
from lib import snake, db, scripts
import os
from sqlalchemy.orm import sessionmaker, aliased

Ranks = ["superkingdom", "phylum", "class", "order",
         "family", "genus", "species"]

def compute_known(filename):
  known = set()
  computed = {}
  if filename and os.path.exists(filename):
    with open(filename) as f:
      for line in f:
        elems = line.rstrip().split("\t")
        known.add(elems[0])
        computed[int(elems[1])] = elems[2:]
  return known, computed

def compute_rank_taxids(summary_file, session, known_results,
                        computed, outfile):
  for line in summary_file:
    if line[0] == "#":
      continue
    elems = line.rstrip().split("\t")
    accession = elems[0]
    if not accession in known_results:
      query_node_id = int(elems[5])
      ids = computed.get(query_node_id, [])
      if not ids:
        try:
          node_id = query_node_id
          prev_node_id = None
          rank_ids = {}
          while prev_node_id != node_id:
            node_q = session.query(NtNode).filter(NtNode.tax_id==node_id)
            node = node_q.one()
            node_id = node.tax_id
            rank_ids[node.rank] = node.tax_id
            prev_node_id = node_id
            node_id = node.parent_tax_id
          for rank in Ranks:
            if rank in rank_ids:
              ids.append(str(rank_ids[rank]))
            else:
              ids.append("None")
        except exc.NoResultFound:
          ids = [elems[1]]
          ids += ["None"]*len(Ranks)
        computed[query_node_id] = ids
      output = [accession] + ids
      outfile.write("\t".join(output)+"\n")
      outfile.flush()

def open_outfile(args):
  if args["--update"]:
    return open(args["--update"], "a")
  else:
    return sys.stdout

def close_outfile(outfile, args):
  if args["--update"]:
    outfile.close()

def main(args):
  known_results, computed = compute_known(args["--update"])
  outfile = open_outfile(args)
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"])
  Session = sessionmaker(bind=engine)
  session = Session()
  compute_rank_taxids(args["<summary>"], session, known_results, computed,
      outfile)
  close_outfile(outfile, args)

def validated(args):
  return scripts.validate(args, db.args_schema,
                          {"--update": Or(None, len),
                           "<summary>": And(str, Use(open))})

if "snakemake" in globals():
  args = snake.args(snakemake, db.snake_args, input=["<summary>"],
                    params=[("--update", "prev")])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(db_args = db.args_doc, common=scripts.args_doc,
    db_args_usage = db.args_usage), version="0.1")
  main(validated(args))
