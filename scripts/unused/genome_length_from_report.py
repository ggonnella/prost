#!/usr/bin/env python3
"""
Download an assembly report and store the genome length into the single_int_stats table.

Usage:
  genomelen_from_report.py <dbuser> <dbpass> <dbname> <dbsocket> <accession>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  accession: NCBI Refseq accession
"""

from sqlalchemy import create_engine
from docopt import docopt
from schema import Schema, And, Use, Or, Optional
import os
import sh
from sqlalchemy.orm import sessionmaker
from dbschema.ncbi_assembly_summary import NcbiAssemblySummary
from dbschema.single_int_stat import SingleIntStat

def main(args):
  connstr = "".join(["mysql+mysqldb://", args["<dbuser>"],
                     ":", args["<dbpass>"], "@localhost/",
                     args["<dbname>"], "?unix_socket=",
                     args["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  Session = sessionmaker(bind=engine)
  session = Session()
  ftp_path = session.query(NcbiAssemblySummary).\
      filter(NcbiAssemblySummary.accession==args["<accession>"]).one().ftp_path
  fn = ftp_path.split("/")[-1]+"_assembly_stats.txt"
  rfn = os.path.join(ftp_path, fn)
  sh.wget(rfn,O=fn)
  glen = int(sh.grep("(?<=^all\\tall\\tall\\tall\\ttotal-length\\t)\d+",fn,
    P=True, o=True))
  os.remove(fn)
  session.add(SingleIntStat(accession=args["<accession>"],
    statname="genome_length", value=glen, source=rfn))
  session.commit()

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<accession>": str,
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<dbuser>": snakemake.config["dbuser"],
    "<dbpass>": snakemake.config["dbpass"],
    "<dbname>": snakemake.config["dbname"],
    "<dbsocket>": snakemake.input.socket,
    "<accession>": snakemake.params.accession}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
