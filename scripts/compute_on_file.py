#!/usr/bin/env python3
"""
Perform computations on a file and output the results as table.

Usage:
  compute_on_file.py [options] <analyzer> <datafile>

Arguments:
  analyzer: py module providing analyze(datafile) => (results, logs)
            analyzes the data file and returns:
            <results>: the analysis results, dict str => anything
                       keys are the measurement name (e.g. "number_of_genes"),
                       values are the measure values (e.g. 12)
            <logs>:    messages reporting unexpected data, errors, etc;
                       keys are categories (e.g. "unexpected_features"),
                       values are lists of strings (e.g. GFF lines)
  datafile: the file to analyze

Options:
  --out, -o FILENAME   output file (default: standard output)
  --log, -l FILENAME   log additional info to the given file
  --pfx, -p INFO       prefix each output line with <INFO> + <TAB>
  --verbose, -v        be verbose
  --version, -V        show script version
  --help, -h           show this help message
"""

from docopt import docopt
from schema import Schema, Use, Optional, Or
import os
from lib import snake, mod, valid

def compute(analyzer, datafile, outfile, logfile, pfx, verbose):
  counters, info = analyzer.analyze(datafile)
  pfx = pfx + "\t" if pfx else ""
  for k, v in counters.items():
    outfile.write(f"{pfx}{k}\t{v}\n")
  if logfile:
    for k, v in info.items():
      for e in v:
        logfile.write(f"{pfx}{k}\t{e}\n")

def main(args):
  analyzer = mod.importer(args["<analyzer>"], args["--verbose"])
  compute(analyzer, args["<datafile>"], args["--out"],
          args["--log"], args["--pfx"], args["--verbose"])

def validated(args):
  schema = Schema({"<datafile>": os.path.exists,
                   "<analyzer>": os.path.exists,
                   "--verbose": Or(None, bool),
                   "--out": valid.outfile_or_stdout,
                   "--log": valid.outfile_or_none,
                   "--pfx": Or(None, str),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        input=["<analyzer>", "<datafile>"],
        output=["--log", "--out"],
        params=["--verbose", "--pfx"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
