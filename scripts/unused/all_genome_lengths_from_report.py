#!/usr/bin/env python3
"""
Download an assembly report and store the genome length into the single_int_stats table.

Usage:
  genomelen_from_report.py <dbuser> <dbpass> <dbname> <dbsocket>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
"""

from sqlalchemy import create_engine
from docopt import docopt
from schema import Schema, And, Use, Or, Optional
import os
import sh
from sqlalchemy.orm import sessionmaker
from dbschema.ncbi_assembly_summary import NcbiAssemblySummary
from dbschema.single_int_stat import SingleIntStat
import tqdm

def main(args):
  connstr = "".join(["mysql+mysqldb://", args["<dbuser>"],
                     ":", args["<dbpass>"], "@localhost/",
                     args["<dbname>"], "?unix_socket=",
                     args["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  Session = sessionmaker(bind=engine)
  session = Session()
  ftp_paths = session.query(NcbiAssemblySummary.accession,
                            NcbiAssemblySummary.ftp_path).\
      filter(NcbiAssemblySummary.seqdb=="refseq").\
      filter(NcbiAssemblySummary.assembly_level=="Complete Genome").\
      filter(NcbiAssemblySummary.ftp_path!="na").\
      filter(NcbiAssemblySummary.version_status=="latest").all()
  for accession, ftp_path in tqdm.tqdm(ftp_paths):
    fn = ftp_path.split("/")[-1]+"_assembly_stats.txt"
    rfn = os.path.join(ftp_path, fn)
    sh.wget(rfn,O=fn)
    glen = int(sh.grep("(?<=^all\\tall\\tall\\tall\\ttotal-length\\t)\d+",fn,
      P=True, o=True))
    os.remove(fn)
    session.add(SingleIntStat(accession=accession,
      statname="genome_length", value=glen, source=rfn))
  session.commit()

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<dbuser>": snakemake.config["dbuser"],
    "<dbpass>": snakemake.config["dbpass"],
    "<dbname>": snakemake.config["dbname"],
    "<dbsocket>": snakemake.input.socket}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
