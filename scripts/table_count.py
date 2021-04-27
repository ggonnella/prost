#!/usr/bin/env python3
"""
Count lines in a tabular file where given sets of conditions are met.

Usage:
  table_count.py [options] <module> [<table>]

Arguments:
  table:  tabular file (standard input if none provided)
          either a text file or gzip-compressed
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
  --check              run module check() function on each row
{common}
"""

from docopt import docopt
from schema import Or
from collections import defaultdict
import os
from lib import snake, valid, tables, mod, scripts

def main(args):
  counts = defaultdict(int)
  m = mod.py(args["<module>"], args["--verbose"])
  for row in tables.get_dict_reader(args, args["<table>"]):
    for k, v in m.counters.items():
      if v(row): counts[k] += 1
    if args["--check"]:
      m.check(row)
  for k in m.counters.keys():
    print(f"{k}{args['--delimiter']}{counts[k]}")

def validated(args):
  return scripts.validate(args, {"<table>": valid.maygz_or_stdin,
                   "<module>": os.path.exists,
                   "--verbose": Or(None, bool),
                   "--comments": valid.comments,
                   "--delimiter": valid.delimiter,
                   "--check": Or(None, bool)})

if "snakemake" in globals():
  args = snake.args(snakemake,
        input=["<table>", "<module>"],
        params=["--verbose", "--fields", "--comments", "--delimiter", "--check"]
      )
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(common=scripts.args_doc), version=0.1)
  main(validated(args))
