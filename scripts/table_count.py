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
import sys
from collections import defaultdict
import importlib
from pathlib import Path
import os
import csv

def import_mod(filename, verbose=False):
  modulename = Path(filename).stem
  spec = importlib.util.spec_from_file_location(modulename, filename)
  m = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(m)
  if verbose:
    sys.stderr.write(f"# module {modulename} imported from file {filename}\n")
  return m

def decomment(f, pfx, verbose=False):
  for line in f:
    if not pfx or not line.startswith(pfx):
      yield line
    elif verbose:
      sys.stderr.write("# skipped internal comment line: "+line)

def get_fieldnames(fieldslist, delimiter, commentspfx, verbose, f):
  if fieldslist:
    fieldnames=fieldslist.split(",")
  else:
    line = f.readline()
    while delimiter not in line:
      if verbose:
        sys.stderr.write("# skipped initial comment line: "+line)
      line = f.readline()
    if line.startswith(commentspfx):
      line = line[len(commentspfx):]
    fieldnames = [x.strip() for x in line.split(delimiter)]
  if verbose:
    sys.stderr.write("# table field names: "+", ".join(fieldnames)+"\n")
  return fieldnames

def main(args):
  counts = defaultdict(int)
  counters = import_mod(args["<module>"], args["--verbose"]).counters
  fieldnames = get_fieldnames(args["--fields"], args["--delimiter"],
                              args["--comments"], args["--verbose"],
                              args["<table>"])
  dictreader = csv.DictReader(decomment(args["<table>"], args["--comments"],
                                        args["--verbose"]),
                            delimiter=args["--delimiter"],
                            fieldnames=fieldnames,
                            quoting=csv.QUOTE_NONE)
  for row in dictreader:
    for k, v in counters.items():
      if v(row): counts[k] += 1
  for k in counters.keys():
    print(f"{k}{args['--delimiter']}{counts[k]}")

def validated(args):
  schema = Schema({"<table>": Or(And(None, Use(lambda v: sys.stdin)), Use(open)),
                   "<module>": os.path.exists,
                   "--verbose": Or(None, bool),
                   "--comments": Or(And(None, Use(lambda v: "#")), str),
                   "--delimiter": Or(And(None, Use(lambda v: "\t")), And(str, len)),
                   Optional(str): object})
  return schema.validate(args)

def args_from_snakemake(args, attr, *names):
  for name in names:
    if name.startswith("--"):
      sname = name[2:]
    else:
      assert(name[0]=="<")
      assert(name[-1]==">")
      sname = name[1:-1]
    sname = sname.replace("-","_")
    args[name] = getattr(snakemake, attr).get(sname)

if "snakemake" in globals():
  args = {}
  args_from_snakemake(args, "input", "<table>", "<module>")
  args_from_snakemake(args, "params", "--verbose", "--fields",
                            "--comments", "--delimiter")
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
