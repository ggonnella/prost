#!/usr/bin/env python3
"""
Compute and output to TSV a stat for all assemblies.

Usage:
  genomestat_for_assembly_files.py [options] <module> <globpattern>

Arguments:
  module:       python module for the stat computation with a function
                value(fn), computing the stat, given a filename
  globpattern:  glob pattern of the directories containing the sequence files

Options
  --skip, -s F     skip accessions in file; this can contain one accession per
                   line, or be a tsv with accession in the first column
  --update, -u F   append output to file (default: output to stdout)
                   (this can be the same file as --skip or another)
{common}
"""

from lib import scripts, snake
from docopt import docopt
from schema import Or
import os
import tqdm
from glob import glob
import sys
import importlib

def main(args):
  spec = importlib.util.spec_from_file_location("stat_comp", args["<module>"])
  stat_comp = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(stat_comp)
  skip = set()
  if args["--skip"] and os.path.exists(args["--skip"]):
    with open(args["--skip"]) as f:
      for line in f:
        skip.add(line.rstrip().split("\t")[0])
  outfile = open(args["--update"], "a") if args["--update"] else sys.stdout
  files = glob(os.path.join(args["<globpattern>"], "*_genomic.fna.gz"))
  for fn in tqdm.tqdm(files):
    accession="_".join(fn.split("/")[-1].split("_")[:2])
    if accession in skip:
      skip.remove(accession)
    else:
      value = stat_comp.value(fn)
      outfile.write("\t".join([accession, value]) + "\n")
  if args["--update"]: outfile.close()

def validated(args):
  return scripts.validate(args, {"<globpattern>": str},
      {"--skip": Or(None, str), "--update": Or(None, str)})

if "snakemake" in globals():
  args = snake.args(snakemake, input=["<module>"], params=["<globpattern>",
    "--skip", "--update"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(common=scripts.args_doc), version=0.1)
  main(validated(args))
