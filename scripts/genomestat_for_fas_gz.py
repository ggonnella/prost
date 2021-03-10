#!/usr/bin/env python3
"""
Compute a stat for a single fas.gz file.

Usage:
  genomestat_for_fas_gz.py [options] <module> <genome>

Arguments:
  module:       python module for the stat computation with a function
                value(fn), computing the stat, given a filename
  genome:       name of the gzipped fasta file

Options
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""

from docopt import docopt
from schema import Schema, Or, Optional
import importlib

def main(args):
  spec = importlib.util.spec_from_file_location("stat_comp", args["<module>"])
  stat_comp = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(stat_comp)
  print(stat_comp.value(args["<genome>"]))

def validated(args):
  schema = Schema({Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<module>": snakemake.input.module,
    "<genome>": snakemake.input.genome,
    }
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
