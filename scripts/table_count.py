#!/usr/bin/env python3
"""
Count lines in a tabular file where given sets of conditions are met.

Usage:
  feature_table_counts.py [options] <module> [<file>]

Arguments:
  file:   tabular file (standard input if none provided)
  module: python module, providing a "count_functions" dictionary
          where the keys are strings and the values are
          callables taking a single argument (line fields)

Output:
  tabular file with two columns (key, count)
  with a row for each key in "count_functions"

Options:
  --delimiter, -d D  field delimiter, default "<TAB>"
  --comments, -c X   prefix of comment lines, default "#"
                     use empty string to keep all lines
  --verbose, -v      be verbose
  --version, -V      show script version
  --help, -h         show this help message
"""

from docopt import docopt
from schema import Schema, Use, Optional, Or, And
import sys
from collections import defaultdict
import importlib
from pathlib import Path
import os

def import_mod(filename):
  modulename = Path(filename).stem
  spec = importlib.util.spec_from_file_location(modulename, filename)
  m = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(m)
  return m

def main(args):
  counts = defaultdict(int)
  count_functions = import_mod(args["<module>"]).count_functions
  for line in args["<file>"]:
    if args["--comments"] and not line.startswith(args["--comments"]):
      elems = line.split(args["--delimiter"])
      for k, v in count_functions.items():
        if v(elems): counts[k] += 1
  for k in count_functions.keys():
    print(f"{k}{args['--delimiter']}{counts[k]}")

def validated(args):
  schema = Schema({"<file>": Or(And(None, Use(lambda v: sys.stdin)), Use(open)),
                   "<module>": os.path.exists,
                   "--comments": Or(And(None, Use(lambda v: "#")), str),
                   "--delimiter": Or(And(None, Use(lambda v: "\t")), And(str, len)),
                   Optional(str): object})
  return schema.validate(args)


if "snakemake" in globals():
  args = {
      "<file>":      snakemake.input.data,
      "<module>":    snakemake.input.module,
      "--comments":  snakemake.params.get("comments"),
      "--delimiter": snakemake.params.get("delimiter"),
    }
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
