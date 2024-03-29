#!/usr/bin/env python3
"""
Count the leaves keys and values from the
output by bacdive_show_structure.py

Usage:
  bacdive_count_leaves.py [options] <file>

Arguments:
  file   file output by bacdive_show_structure.py

Options:
{common}
"""

from lib import scripts
from schema import Use
import snacli

def perc(n_part, n_all):
  return "{:.1f}%".format((n_part / n_all) * 100)

def main(args):
  ids = set()
  for line in args["<file>"]:
    elems = line.rstrip().split("\t")
    ids.add(elems[2])
  n_strains = len(ids)
  leaves = dict()
  args["<file>"].seek(0)
  for line in args["<file>"]:
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

def validated(args):
  return scripts.validate(args, {"<file>": Use(open)})

with snacli.args(docvars={"common": scripts.args_doc},
                 input=[("<file>", "data")], version="0.1") as args:
   if args: main(validated(args))
