#!/usr/bin/env python3
"""
Extract all accessions from the Bacdive strain
information TSV produced by download_strain_info.py.

Usage:
  ./extract_bacdive_accessions.py [options] <tsv>

Arguments:
  tsv   output of download_strain_info.py

Options:
  --verbose, -v   be verbose
  --help, -h      show this help message
  --version, -V   show the script version
"""

from docopt import docopt
import json

def traverse(j, bacdive_id, pfx, last):
  if isinstance(j, dict):
    for k, v in j.items():
      traverse(v, bacdive_id, pfx + "\\" + k, k)
  elif isinstance(j, list):
    for v in j:
      traverse(v, bacdive_id, pfx+"[]", last+"[]")
  elif last == "seq_acc_num":
    print(j + "\t" + bacdive_id)

def main(arguments):
  with open(arguments["<tsv>"]) as f:
    for line in f:
      elems = line.rstrip().split("\t")
      bacdive_id = elems[0]
      details = json.loads(elems[3])
      traverse(details, bacdive_id, "", "")

if __name__ == "__main__":
  arguments = docopt(__doc__, version=0.1)
  main(arguments)
