#!/usr/bin/env python3
"""
Match accessions in two TSV files

Usage:
  match_accessions.py [options] <t1> <m1> <o1> <t2> <m2> <o2>

Arguments:
  t1: first TSV file
  m1: column w. accession to match, in file t1
  o1: column to output from file t1
  t2: second TSV file
  m2: column w. accession to match, in file t2
  o2: column to output from file t2
  column numbers (m1, m2) are 1-based

Options:
  --strip, -s    strip version number from accessions before match
  --reverse, -r  output column o2 before o1
  --multi, -m    column m2 contains multiple accessions, comma-separated
{common}
"""

from schema import Use, And
from lib import scripts
import snacli

def strip_version_number(acc_v):
  return acc_v.split(".")[0]

def main(args):
  t1data = dict()
  for line in args["<t1>"]:
    elems = line.rstrip().split("\t")
    key = elems[args["<m1>"]]
    if args["--strip"]:
      key = strip_version_number(key)
    value = elems[args["<o1>"]]
    t1data[key]=value
  for line in args["<t2>"]:
    elems = line.rstrip().split("\t")
    keys = elems[args["<m2>"]]
    if args["--multi"]:
      keys = keys.split(",")
    else:
      keys = [keys]
    for key in keys:
      if args["--strip"]:
        key = strip_version_number(key)
      if key in t1data:
        value1 = t1data[key]
        value2 = elems[args["<o2>"]]
        if args["--reverse"]:
          print(f"{value2}\t{value1}")
        else:
          print(f"{value1}\t{value2}")

def validated(args):
  colnum = And(Use(lambda i: int(i)-1), lambda i: i>=0)
  return scripts.validate(args, {
    "<t1>": Use(open), "<t2>": Use(open),
    "<m1>": colnum, "<m2>": colnum, "<o1>": colnum, "<o2>": colnum})

with snacli.args(input=["<t1>", "<t2>"],
                 params=["<m1>", "<m2>", "<o1>", "<o2>"],
                 docvars={"common": scripts.args_doc},
                 version="0.1") as args:
  if args: main(validated(args))
