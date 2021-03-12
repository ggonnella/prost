#!/usr/bin/env python3
"""
Perform computations on a GFF file and output the results as table.

Usage:
  gff_compute.py [options] <module> <gff>

Arguments:
  gff:    GFF file
  module: python module, providing a results(db) function
          where db is a gffutils database object;
          the function returns a dict of results
          with string key and values

Options:
  --verbose, -v        be verbose
  --version, -V        show script version
  --help, -h           show this help message
"""

from docopt import docopt
from schema import Schema, Use, Optional, Or
import gffutils
import os
from lib import snake, mod

def main(args):
  m = mod.importer(args["<module>"], args["--verbose"])
  db = gffutils.create_db(args["<gff>"], ":memory:",
                          merge_strategy="create_unique")
  for k,v in m.results(db).items():
    print(f"{k}\t{v}")

def validated(args):
  schema = Schema({"<gff>": os.path.exists,
                   "<module>": os.path.exists,
                   "--verbose": Or(None, bool),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        input=["<gff>", "<module>"],
        params=["--verbose"]
      )
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
