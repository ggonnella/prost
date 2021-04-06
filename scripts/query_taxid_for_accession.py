#!/usr/bin/env python3
"""
Output taxid for the given NCBI assembly accession.

Usage:
  query_taxid_for_accession.py [options] {db_args_usage}
                               <accession>

Arguments:
{db_args}
  accession: NCBI assembly accession

Options:
{common}
"""

from sqlalchemy import create_engine
from dbschema.ncbi_assembly_summary import NcbiAssemblySummary
from docopt import docopt
from lib import db, snake, scripts
from sqlalchemy.orm import sessionmaker

def main(args):
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"])
  Session = sessionmaker(bind=engine)
  session = Session()
  taxid = session.query(NcbiAssemblySummary.taxid).\
      filter(NcbiAssemblySummary.accession==args["<accession>"]).one()
  print(taxid[0])

def validated(args):
  return scripts.validate(args, db.args_schema, {"<accession>": str})

if "snakemake" in globals():
  args = snake.args(snakemake, db.snake_args, params = ["<accession>"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(db_args = db.args_doc, common = scripts.args_doc,
    db_args_usage = db.args_usage), version="0.1")
  main(validated(args))
