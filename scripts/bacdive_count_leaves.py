#!/usr/bin/env python3
"""
Count the leaves keys and values from the
output by bacdive_show_structure.py

Usage:
  bacdive_count_leaves.py [options] <file>

Arguments:
  file   file output by bacdive_show_structure.py

Options:
  --verbose, -v  be verbose
  --version, -V  show script version
  --help, -h     show this help message
"""

from docopt import docopt

def perc(n_part, n_all):
  return "{:.1f}%".format((n_part / n_all) * 100)

def main(arguments):
  ids = set()
  with open(arguments["<file>"]) as f:
    for line in f:
      elems = line.rstrip().split("\t")
      ids.add(elems[2])
  n_strains = len(ids)
  leaves = dict()
  with open(arguments["<file>"]) as f:
    for line in f:
      elems = line.rstrip().split("\t")
      key = elems[1]
      value = elems[3]
      if key not in leaves:
        leaves[key] = dict()
      if value not in leaves[key]:
        leaves[key][value] = 1
      else:
        leaves[key][value] += 1
  for k, k_stats in leaves.items():
    total = 0
    for v, c in k_stats.items():
      total += c
    for v, c in k_stats.items():
      print(f"{k}\t{v}\t{c}\t{perc(c, total)}")
    print(f"#{k}\ttotal\t{total}\t{perc(total, n_strains)}")
    print("")

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(arguments)
