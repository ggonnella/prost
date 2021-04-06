#!/usr/bin/env python3
"""
Remove sequence which were downloaded twice, because of programming errors.

Usage:
  remove_dup_seq.py [options] <globpattern>

Arguments:
  globpattern: glob pattern of the directories containing chunk directories
               (chunk.*) containing the sequence files (*_genomic.fna.gz)
               i.e. "/chunk.*/*_genomic.fna.gz" is appended to it
Options:
{common}
"""
from docopt import docopt
from glob import glob
from lib import scripts, snake
import os
import sys

def path2accession(path):
  return "_".join(path.split("/")[-1].split("_")[:2])

def path2chunknum(path):
  return int(path.split("/")[-2].split(".")[1])

def main(args):
  allpaths = glob(args["<globpattern>"]+"/chunk.*/*_genomic.fna.gz")
  accession2chunk = {}
  for path in allpaths:
    accession = path2accession(path)
    chunknum = path2chunknum(path)
    if accession not in accession2chunk:
      accession2chunk[accession] = {}
    accession2chunk[accession][chunknum] = path
  for accession, chunknum2path in accession2chunk.items():
    if len(chunknum2path) > 1:
      keep = min(chunknum2path.keys())
      sz_keep = os.stat(chunknum2path[keep]).st_size
      for chunknum, path in chunknum2path.items():
        if chunknum != keep:
          sz_f = os.stat(path).st_size
          if sz_f == sz_keep:
            sys.stderr.write(f"removing {path} since it is the duplicate "+
                             f"of {chunknum2path[keep]}\n")
            os.remove(path)
          else:
            sys.stderr.write(f"not removing {path} since it is the duplicate "+
                             f"of {chunknum2path[keep]} but file sizes differ "+\
                             f"({sz_f} != {sz_keep})\n")

def validated(args):
  return scripts.validate(args, {"<globpattern>": str})

if "snakemake" in globals():
  args = snake.args(snakemake, params = ["<globpattern>"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(common=scripts.args_doc), version=0.1)
  main(validated(args))
