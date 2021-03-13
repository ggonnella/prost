#!/usr/bin/env python3
"""
Perform computations on a GFF file and output the results as table.

Usage:
  gff_compute.py [options] <module> <gff>

Arguments:
  gff:    GFF file (uncompressed or gzipped)
  module: python module, which must provide a function:
            analyze(db) => (counters, info)
          <db>: "gffutils" database object
          <counters>: a dict of strings => int
                      where expected features are counted
          <info>: a dict of string => list of strings
                  where unexpected features or errors are stored

Options:
  --out, -o FILE       output file (default: standard output)
  --info, -i INFO      log the info to the given file
  --verbose, -v        be verbose
  --version, -V        show script version
  --help, -h           show this help message
"""

from docopt import docopt
from schema import Schema, Use, Optional, Or
import gffutils
import os
from lib import snake, mod, valid

def main(args):
  m = mod.importer(args["<module>"], args["--verbose"])
  db = gffutils.create_db(args["<gff>"], ":memory:",
                          merge_strategy="create_unique")
  counters, info = m.analyze(db)
  for k, v in counters.items():
    args["--out"].write(f"{k}\t{v}\n")
  if args["--info"]:
    for k, v in info.items():
      for e in v:
        args["--info"].write(f"{k}\t{e}\n")

def validated(args):
  schema = Schema({"<gff>": os.path.exists,
                   "<module>": os.path.exists,
                   "--verbose": Or(None, bool),
                   "--out": valid.outfile_or_stdout,
                   "--info": valid.outfile_or_none,
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        input=["<gff>", "<module>"],
        output=["--info", "--out"],
        params=["--verbose"]
      )
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
