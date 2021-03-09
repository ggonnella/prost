#!/usr/bin/env python3
"""
Update the table of Refseq assembly accessions to Refseq sequence accessions.

Usage:
  rs_seq_accessions_for_assemblies.py [options] <globpattern>

Arguments:
  globpattern: glob pattern to locate the directories contianing the
               Refseq assemblies files

Options:
  --update, -u F            file with previous results, to update
  --verbose, -v             be verbose
  --version, -V             show script version
  --help, -h                show this help message
"""

from docopt import docopt
import sh
import os
import sys
import tqdm
import glob
from schema import Schema, Use, Optional, Or, And

def extract_fastaids(localfn):
  # zcat localfn | grep -P '^>' | cut -f1 -d' ' | cut -c2-
  seqaccs=sh.cut(sh.cut(sh.grep(sh.zcat(localfn,_piped=True),
    "^>",P=True,_piped=True), f="1",d=" ",_piped=True),
    c="2-").split("\n")[:-1]
  return seqaccs

def compute_known(filename):
  result = set()
  if filename and os.path.exists(filename):
    with open(filename) as f:
      for line in f:
        elems = line.rstrip().split("\t")
        result.add(elems[0])
  return result

def compute_targets(globpattern, known_results):
  result = []
  for fn in glob.glob(os.path.join(globpattern, "*_genomic.fna.gz")):
    asmacc = "_".join(os.path.basename(fn).split("_")[:2])
    if asmacc not in known_results:
      result.append((asmacc, fn))
  return result

def compute_accessions(targets, outfile):
  noferrors = 0
  for asmacc, fn in tqdm.tqdm(targets):
    sys.stderr.write(f"Processing {fn}...\n")
    seqaccs = extract_fastaids(fn)
    outfile.write("\t".join([asmacc, ",".join(seqaccs)])+"\n")
    outfile.flush()

def open_outfile(args):
  if args["--update"]:
    return open(args["--update"], "a")
  else:
    return sys.stdout

def close_outfile(outfile, args):
  if args["--update"]:
    outfile.close()

def main(args):
  known_results = compute_known(args["--update"])
  targets = compute_targets(args["<globpattern>"], known_results)
  outfile = open_outfile(args)
  compute_accessions(targets, outfile)
  close_outfile(outfile, args)

def validated(args):
  schema = Schema({
    Optional("--update"): Or(None, os.path.exists),
    "<globpattern>": str,
    Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
      "--update": snakemake.params.get("prev", None),
      "<globpattern>": snakemake.params.globpattern,
    }
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
