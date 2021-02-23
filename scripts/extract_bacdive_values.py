#!/usr/bin/env python3
"""
Extract all values for a given key from the Bacdive strain
information TSV produced by download_strain_info.py.

Usage:
  ./extract_bacdive_values.py [options] <tsv> <key>

Arguments:
  tsv   output of download_strain_info.py
  key   key for which the values shall be extracted

Options:
  --split, -s D   split value by delimiter
  --verbose, -v   be verbose
  --help, -h      show this help message
  --version, -V   show the script version
"""

from docopt import docopt
from schema import Schema, Use, Or
import json

def traverse(j, bacdive_id, pfx, last, key, split):
  if isinstance(j, dict):
    for k, v in j.items():
      traverse(v, bacdive_id, pfx + "\\" + k, k, key, split)
  elif isinstance(j, list):
    for v in j:
      traverse(v, bacdive_id, pfx+"[]", last+"[]", key, split)
  elif last == key:
    if j:
      if split:
        values = [v.strip() for v in j.split(split)]
      else:
        values = [j]
      for jelem in values:
        print(jelem + "\t" + bacdive_id)

def main(arguments):
  for line in arguments["<tsv>"]:
      elems = line.rstrip().split("\t")
      bacdive_id = elems[0]
      details = json.loads(elems[3])
      traverse(details, bacdive_id, "", "", arguments["<key>"],
          arguments["--split"])

def validated(arguments):
  schema = Schema({
    "<tsv>": Use(open), "<key>": str,
    "--split": Or(None, str)
    }, ignore_extra_keys=True)
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version=0.1)
  main(validated(arguments))
