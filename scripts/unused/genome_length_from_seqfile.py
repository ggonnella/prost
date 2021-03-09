#!/usr/bin/env python3
"""
Download an assembly report and store the genome length into the single_int_stats table.

Usage:
  genomelen_from_report.py <dbuser> <dbpass> <dbname> <dbsocket> <globpattern> <accession>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  globpattern: glob pattern of the directories containing the sequence files
  accession: NCBI Refseq accession
"""

from sqlalchemy import create_engine
from docopt import docopt
from schema import Schema, And, Use, Or, Optional
import os
import sh
from sqlalchemy.orm import sessionmaker
from glob import glob
from dbschema.single_int_stat import SingleIntStat

def main(args):
  connstr = "".join(["mysql+mysqldb://", args["<dbuser>"],
                     ":", args["<dbpass>"], "@localhost/",
                     args["<dbname>"], "?unix_socket=",
                     args["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  Session = sessionmaker(bind=engine)
  session = Session()
  files = glob(os.path.join(args["<globpattern>"],
               args["<accession>"]+"*_genomic.fna.gz"))
  assert(len(files) == 1)
  fn = files[0]
  glen = sh.wc(sh.tr(sh.grep(sh.zcat(fn, _piped=True), ">", v=True, _piped=True),
    d="[:space:]", _piped=True), c=True)
  session.add(SingleIntStat(accession=args["<accession>"],
    statname="genome_length", value=glen, source="computed from fna.gz"))
  session.commit()

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<accession>": str,
                   "<globpattern>": str,
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<dbuser>": snakemake.config["dbuser"],
    "<dbpass>": snakemake.config["dbpass"],
    "<dbname>": snakemake.config["dbname"],
    "<dbsocket>": snakemake.input.socket,
    "<globpattern>": snakemake.params.globepattern,
    "<accession>": snakemake.params.accession}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
