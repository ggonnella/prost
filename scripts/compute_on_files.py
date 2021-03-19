#!/usr/bin/env python3
"""
Perform computations on multiple files and output the results as table.

Usage:
  compute_on_file.py [options] <analyzer> <globpattern>

Arguments:
  analyzer: py module providing analyze(datafile, **kwargs) => (results, logs)
            analyzes the data file and returns:
            <results>: the analysis results, dict {str: dict {str: any}}
                       keys of external dict: attr_class (e.g. "seq_stats")
                       keys of internal dict: attr_instance (e.g. "genome_size")
                       values: attribute value (e.g. 12)
            <logs>:    messages reporting unexpected data, errors, etc;
                       keys are categories (e.g. "unexpected_features"),
                       values are lists of strings (e.g. GFF lines)
            kwargs are optionally additional parameters for the computation
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
  --report, -r FN   computation report file (default: stderr)
  --user U          user_id for the report (default: getpass.getuser())
  --system S        system_id for the report (default: socket.gethostname())
  --reason R        reason field for the report (default: None)
  --params FNAME    YAML file with additional parameters (default: None)
  --verbose, -v     be verbose
  --version, -V     show script version
  --help, -h        show this help message
"""

from docopt import docopt
from schema import Schema, Use, Optional, Or, And
import os
import sys
from glob import glob
from lib import snake, mod, valid
import tqdm
import socket
import getpass
import yaml
from datetime import datetime
import textwrap
import uuid

def compute(analyzer, datafile, params, outfile, logfile, pfx, verbose):
  if analyzer.__nim__:
    counters, info = analyzer.analyze(datafile, params)
  else:
    counters, info = analyzer.analyze(datafile, **params)
  pfx = pfx + "\t" if pfx else ""
  for k1, k2_v in counters.items():
    for k2, v in k2_v.items():
      outfile.write(f"{pfx}{k1}\t{k2}\t{v}\n")
  if logfile and info:
    for k, v in info.items():
      for e in v:
        logfile.write(f"{pfx}{k}\t{e}\n")

class Report():
  def __init__(self, analyzer, user, system, reason, params):
    self.data = {}
    plugin_data = yaml.safe_load(analyzer.analyze.__doc__)
    self.data["plugin_id"] = plugin_data["id"]
    self.data["plugin_version"] = plugin_data["version"]
    self.data["system_id"] = system
    self.data["user_id"] = user
    if reason:
      self.data["reason"] = reason
    if params:
      self.data["parameters"] = yaml.dump(params)
    self.data["uuid"] = uuid.uuid4().bytes
    self.data["time_start"] = str(datetime.now())

  def finalize(self, rfile, n_processed, err=None, datafile=None):
    self.data["time_end"] = str(datetime.now())
    self.data["n_units"] = n_processed
    if err:
      self.data["status"] = "aborted" if n_processed == 0 else "partial"
      remark = {}
      remark["error_datafile"] = datafile
      remark["error_class"] = err.__class__.__name__
      remark["error_message"] = str(err)
      self.data["remarks"] = yaml.dump(remark)
    else:
      self.data["status"] = "completed"
    yaml.dump(self.data, rfile)
    rfile.flush()

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
  analyzer = mod.py_or_nim(args["<analyzer>"], args["--verbose"])
  report = Report(analyzer, args["--user"], args["--system"], args["--reason"],
                  args["--params"])
  outfile = open(args["--out"], "a") if args["--out"] else sys.stdout
  logfile = open(args["--log"], "a") if args["--log"] else sys.stderr
  files = glob(args["<globpattern>"])
  n_processed = 0
  for fn in tqdm.tqdm(files):
    identifier = id_from_filename(fn)
    if identifier in skip:
      skip.remove(identifier)
    else:
      try:
        compute(analyzer, fn, args["--params"],
                outfile, logfile, identifier, args["--verbose"])
        n_processed += 1
      except Exception as err:
        logfile.flush()
        outfile.flush()
        report.finalize(args["--report"], n_processed, err, fn)
        raise(err)
  report.finalize(args["--report"], n_processed)
  if args["--out"]: outfile.close()
  if args["--log"]: logfile.close()

def validated(args):
  schema = Schema({"--fnparser": Or(None, os.path.exists),
                   "<globpattern>": str,
                   "<analyzer>": os.path.exists,
                   "--verbose": Or(None, bool),
                   "--out": Or(None, str),
                   "--log": Or(None, str),
                   "--report": valid.outfile_or_stderr,
                   "--user": Or(str, And(None,
                                 Use(lambda n: getpass.getuser()))),
                   "--system": Or(str, And(None,
                                   Use(lambda n: socket.gethostname()))),
                   "--params": Or(And(None, Use(lambda n: {})),
                                  And(str, Use(open),
                                    Use(lambda fn: yaml.safe_load(fn)))),
                   "--skip": Or(None, os.path.exists),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        input=["<analyzer>", "--fnparser"],
        log=["--out", "--log"],
        params=["--verbose", "<globpattern>",
                "--user", "--system", "--reason"])
  if args["--out"] and os.path.exists(args["--out"]):
    args["--skip"] = args["--out"]
  else:
    args["--skip"] = None
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
