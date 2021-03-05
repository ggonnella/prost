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
  --verbose, -v  be verbose
  --help, -h     show this help message
  --version, -V  show the script version
  --multi, -m    column m2 contains multiple accessions, comma-separated
"""

from docopt import docopt
from schema import Schema, Use, And, Optional
import sys

def strip_version_number(acc_v):
  return acc_v.split(".")[0]

def main(arguments):
  t1data = dict()
  for line in arguments["<t1>"]:
    elems = line.rstrip().split("\t")
    key = elems[arguments["<m1>"]]
    if arguments["--strip"]:
      key = strip_version_number(key)
    value = elems[arguments["<o1>"]]
    t1data[key]=value
  for line in arguments["<t2>"]:
    elems = line.rstrip().split("\t")
    keys = elems[arguments["<m2>"]]
    if arguments["--multi"]:
      keys = keys.split(",")
    else:
      keys = [keys]
    for key in keys:
      if arguments["--strip"]:
        key = strip_version_number(key)
      if key in t1data:
        value1 = t1data[key]
        value2 = elems[arguments["<o2>"]]
        if arguments["--reverse"]:
          print(f"{value2}\t{value1}")
        else:
          print(f"{value1}\t{value2}")

def validated(arguments):
  colnum = And(Use(lambda i: int(i)-1), lambda i: i>=0)
  schema = Schema({
    "<t1>": Use(open), "<t2>": Use(open),
    "<m1>": colnum, "<m2>": colnum, "<o1>": colnum, "<o2>": colnum,
    Optional(str): object})
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version=0.1)
  main(validated(arguments))
