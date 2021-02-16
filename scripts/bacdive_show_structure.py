#!/usr/bin/env python3
"""
Show the structure of the data of Bacdive JSON data
for strains.

Usage:
  ./bacdive_show_structure.py [options] <bacdivefile>

Arguments:
  bacdivefile  TSV file containing strain data

Options:
  --help, -h     show this help message
  --version, -V  show the script version
  --verbose, -v  be verbose
"""

from docopt import docopt
import json

def print_leaves(j, bacdive_id, pfx, last):
  if isinstance(j, dict):
    for k, v in j.items():
      print_leaves(v, bacdive_id, pfx + "\\" + k, k)
  elif isinstance(j, list):
    for v in j:
      print_leaves(v, bacdive_id, pfx+"[]", last+"[]")
  else:
    print(pfx + "\t" + last + "\t" + bacdive_id + "\t" + json.dumps(j))

def main(arguments):
  with open(arguments["<bacdivefile>"]) as f:
    for line in f:
      elems = line.rstrip().split("\t")
      bacdive_id = elems[0]
      details = json.loads(elems[3])
      print_leaves(details, bacdive_id, "", "")

if __name__ == "__main__":
  arguments = docopt(__doc__, version=0.1)
  main(arguments)
