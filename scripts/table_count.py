#!/usr/bin/env python3
"""
Count lines in a tabular file where given sets of conditions are met.

Usage:
  table_count.py [options] <module> [<table>]

Arguments:
  table:  tabular file (standard input if none provided)
  module: python module, providing a dictionary "counters"
          keys:   string
          values: callable(row_dictionary) => bool

Output:
  tabular file with two columns (key, count)
  with a row for each key in "counters"

Options:
  --delimiter, -d D    field delimiter, default "<TAB>"
  --comments, -c X     prefix of comment lines, default "#"
                       use empty string to keep all lines
  --fields, -f NAMES   comma-separated list of field names
                       default: from first line containing the delimiter
                                (after removing the comment pfx, if any)
  --verbose, -v        be verbose
  --version, -V        show script version
  --help, -h           show this help message
"""

from docopt import docopt
from schema import Schema, Use, Optional, Or, And
from collections import defaultdict
import os
from lib import snake, valid, tables, mod

def main(args):
  counts = defaultdict(int)
  counters = mod.importer(args["<module>"], args["--verbose"]).counters
  for row in tables.get_dict_reader(args, args["<table>"]):
    for k, v in counters.items():
      if v(row): counts[k] += 1
  for k in counters.keys():
    print(f"{k}{args['--delimiter']}{counts[k]}")

def validated(args):
  schema = Schema({"<table>": valid.open_file_or_stdin,
                   "<module>": os.path.exists,
                   "--verbose": Or(None, bool),
                   "--comments": valid.comments,
                   "--delimiter": valid.delimiter,
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        input=["<table>", "<module>"],
        params=["--verbose", "--fields", "--comments", "--delimiter"]
      )
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
