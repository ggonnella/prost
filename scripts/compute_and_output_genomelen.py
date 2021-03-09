#!/usr/bin/env python3
"""
Compute and output to TSV the genome length for all assemblies.

Usage:
  compute_and_output_genomelen.py [options] <globpattern>

Arguments:
  globpattern:  glob pattern of the directories containing the sequence files

Options
  --skip, -s F     skip accessions in file; this can contain one accession per
                   line, or be a tsv with accession in the first column
  --update, -u F   append output to file (default: output to stdout)
                   (this can be the same file as --skip or another)
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""

from docopt import docopt
from schema import Schema, Or, Optional
import os
import sh
import tqdm
from glob import glob
import sys

def main(args):
  skip = set()
  if args["--skip"] and os.path.exists(args["--skip"]):
    with open(args["--skip"]) as f:
      for line in f:
        skip.add(line.rstrip().split("\t")[0])
  outfile = open(args["--update"], "a") if args["--update"] else sys.stdout
  files = glob(os.path.join(args["<globpattern>"], "*_genomic.fna.gz"))
  statname = "genome_length"
  source = "computed from genomic.fas.gz; obtained by FTP from NCBI Refseq"
  for fn in tqdm.tqdm(files):
    accession="_".join(fn.split("/")[-1].split("_")[:2])
    if accession in skip:
      skip.remove(accession)
    else:
      value = sh.wc(sh.tr(sh.grep(sh.zcat(fn, _piped=True), ">", v=True,
        _piped=True), d="[:space:]", _piped=True), c=True)
      outfile.write("\t".join([accession, statname,
                               str(value).rstrip(), source]) + "\n")
  if args["--update"]: outfile.close()

def validated(args):
  schema = Schema({Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = {
    "<globpattern>": snakemake.params.globpattern,
    "--skip": snakemake.params.get("skip", None),
    "--update": snakemake.params.get("update", None)}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
