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
{common}
"""

from docopt import docopt
import importlib
import os
from lib import snake, scripts

def compute_value(modulename, genomefn):
  spec = importlib.util.spec_from_file_location("stat_comp", modulename)
  stat_comp = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(stat_comp)
  return stat_comp.value(genomefn)

def main(args):
  print(compute_value(args["<module>"], args["<genome>"]))

def validated(args):
  return scripts.validate(args, {"module": os.path.exists,
                                 "genome": os.path.exists})

if "snakemake" in globals():
  args = snake.args(snakemake, input=["<module>", "<genome>"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(common=scripts.args_doc), version=0.1)
  main(validated(args))
