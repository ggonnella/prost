#!/usr/bin/env python3
"""
Output taxid for the given NCBI assembly accession.

Usage:
  query_taxid_for_accession.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <accession>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  accession: NCBI assembly accession

Options:
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""

from sqlalchemy import create_engine
from dbschema.ncbi_assembly_summary import NcbiAssemblySummary
from docopt import docopt
from schema import Schema, And, Use, Or, Optional
import os
from sqlalchemy.orm import sessionmaker, aliased

def main(args):
  connstr = "".join(["mysql+mysqldb://", args["<dbuser>"],
                     ":", args["<dbpass>"], "@localhost/",
                     args["<dbname>"], "?unix_socket=",
                     args["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  Session = sessionmaker(bind=engine)
  session = Session()
  taxid = session.query(NcbiAssemblySummary.taxid).\
      filter(NcbiAssemblySummary.accession==args["<accession>"]).one()
  print(taxid[0])

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
    "<accession>": snakemake.params.accession,
    }
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
