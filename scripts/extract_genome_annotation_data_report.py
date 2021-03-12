#!/usr/bin/env python3
"""
Extract the genome annotation data reports from a Refseq Genbank file.

Usage:
  extract_genome_annotation_data_report.py [options] [<file>]

Arguments:
  file: Refseq Genbank file (stdin if none provided)

Options:
  --one, -1        output only first report found
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""

from docopt import docopt
from schema import Schema, Use, Optional, Or, And
from lib import valid

def main(args):
  state = 0
  leading_spaces = 0
  for line in args["<file>"]:
    line = line.rstrip()
    if line.endswith("##Genome-Annotation-Data-START##"):
      leading_spaces = len(line)-len(line.lstrip())
      state = 1
    elif line.endswith("##Genome-Annotation-Data-END##"):
      state = 2
    elif state == 2:
      state = 0
    if state: print(line[leading_spaces:])
    if state == 2 and args["--one"]:
      exit(0)

def validated(args):
  schema = Schema({"<file>": valid.open_maygz_or_stdin,
                   "--one": Or(None, bool),
                   Optional(str): object})
  return schema.validate(args)


if "snakemake" in globals():
  args = { "<file>": snakemake.input.gb,
           "--one": snakemake.params.one}
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
