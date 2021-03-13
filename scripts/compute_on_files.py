#!/usr/bin/env python3
"""
Perform computations on multiple files and output the results as table.

Usage:
  compute_on_file.py [options] <analyzer> <globpattern>

Arguments:
  analyzer: py module providing analyze(datafile) => (results, logs)
            analyzes the data file and returns:
            <results>: the analysis results, dict str => anything
                       keys are the measurement name (e.g. "number_of_genes"),
                       values are the measure values (e.g. 12)
            <logs>:    messages reporting unexpected data, errors, etc;
                       keys are categories (e.g. "unexpected_features"),
                       values are lists of strings (e.g. GFF lines)
  globpattern: the globpattern to find the files

Options:
  --fnparser FNAME  py module providing a id_from_filename(filename) function,
                    which computes an identifier (e.g. accession);
                    the identifier is written in the first column of the output
                    (default: the identifier is the filename)
  --skip, -s FNAME  skip computations for which the identifier
                    (filename or result of fnparser) is contained in this file
                    (one ID per line, or tsv with IDs in first column)
  --out, -o FNAME   output results to file (default: stdout);
                    if the file exists, the output is appended;
                    (this can be the same file used for --skip)
  --log, -l FNAME   write logs to the given file (default: stderr);
                    if the file exists, the output is appended
  --verbose, -v     be verbose
  --version, -V     show script version
  --help, -h        show this help message
"""

from docopt import docopt
from schema import Schema, Use, Optional, Or
import os
import sys
from glob import glob
from lib import snake, mod, valid
import tqdm
import compute_on_file

def id_from_filename(filename):
  return "_".join(filename.split("/")[-1].split("_")[:2])

def main(args):
  skip = set()
  if args["--skip"]:
    with open(args["--skip"]) as f:
      for line in f: skip.add(line.rstrip().split("\t")[0])
  if args["--fnparser"]:
    fnparser = mod.importer(args["--fnparser"], args["--verbose"])
    id_from_filename = fnparser.id_from_filename
  else:
    id_from_filename = lambda x: x
  analyzer = mod.importer(args["<analyzer>"], args["--verbose"])
  outfile = open(args["--out"], "a") if args["--out"] else sys.stdout
  logfile = open(args["--log"], "a") if args["--log"] else sys.stderr
  files = glob(args["<globpattern>"])
  for fn in tqdm.tqdm(files):
    identifier = id_from_filename(fn)
    if identifier in skip:
      skip.remove(identifier)
    else:
      compute_on_file.compute(analyzer, fn, outfile,
                              logfile, identifier, args["--verbose"])
  if args["--out"]: outfile.close()
  if args["--log"]: logfile.close()

def validated(args):
  schema = Schema({"--fnparser": Or(None, os.path.exists),
                   "<globpattern>": str,
                   "<analyzer>": os.path.exists,
                   "--verbose": Or(None, bool),
                   "--out": Or(None, str),
                   "--log": Or(None, str),
                   "--skip": Or(None, os.path.exists),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        input=["<analyzer>", "--fnparser"],
        log=["--out", "--log"],
        params=["--verbose", "--pfx", "<globpattern>"])
  if args["--out"] and os.path.exists(args["--out"]):
    args["--skip"] = args["--out"]
  else:
    args["--skip"] = None
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
