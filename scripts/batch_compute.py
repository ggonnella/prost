#!/usr/bin/env python3
"""
Perform computations on multiple files, using the compute function of the
specified Nim or Python plugin module (see plugins/README.md for the plugins
specification).

Usage:
  batch_compute.py <plugin> ids <idsfile> [<col>] [options]
  batch_compute.py <plugin> files <globpattern> [options]

The compute function is applied to each input ID or filename (see below).
Unless --skip is used and the input ID (or ID computed from the filename)
is present in the skip file.

Input:

  ids <idsfile> [<col>]

    (1a) IDs in a file, one per line => use "ids <idsfile>"
    (1b) IDs in a TSV file => use "ids <idsfile> <col>"; 1-based col number

  => to edit the ids (passed to compute() and used by --skip and
     in the results) use --idsproc

  files <globpattern>

    (2) Names of files matching a glob pattern => "files <globpattern>"

  => to compute the IDs (used by --skip and in the results), use --idsproc,
     otherwise the filenames are used as IDs;
     (note: the input to compute() are in this case the filenames,
      thus not affected by --idsproc)

Output:
- computation results: TSV file, where the first column contains the IDs,
  the following columns are the results in the order specified by the plugin
- logs: TSV file, first column are the unit IDs, followed by information
  messages returned by the plugin (each computation can generate zero, one
  or multiple lines)
- computation report: YAML file, computation time, status, plugin name,
  version, parameters, username, hostname, etc.

Options:
  --idsproc FNAME   Python or Nim module, providing compute_id(str)->str;
                    allows to edit the IDs/filenames used for (1) results;
                    (2) --skip option; (3) in "ids" mode only: input to the
                    compute() function
  --skip, -s FNAME  skip computations for which the ID is contained in this
                    file (one ID per line, or TSV with IDs in first column)
  --out, -o FNAME   output results to file (default: stdout);
                    if the file exists, the output is appended;
                    (this can be the same file used for --skip)
  --log, -l FNAME   write logs to the given file (default: stderr);
                    if the file exists, the output is appended
{report_opts}
{common}
"""

from docopt import docopt
from schema import Or
import os
import sys
from glob import glob
from lib import snake, mod, valid, reports, scripts
import tqdm

def compute(plugin, state, input_id, output_id, params,
            outfile, logfile, report, n_processed, verbose):
  if state:
    params["state"] = state
  try:
    results, logs = plugin.compute(input_id, **params)
  except Exception as err:
    report.finalize(n_processed, err, input_id)
    outfile.flush()
    logfile.flush()
    raise(err)
  results = "\t".join([str(r) for r in results])
  outfile.write(f"{output_id}\t{results}\n")
  if logfile and logs:
    for msg in logs:
      logfile.write(f"{output_id}\t{msg}\n")

def input_units(args):
  if args["<globpattern>"]:
    return glob(args["<globpattern>"])
  else:
    with open(args["<idsfile>"]) as f:
      return [line.rstrip().split("\t")[args["<col>"]-1] for line in f]

def compute_skip_set(skip_arg, verbose):
  skip = set()
  if skip_arg:
    if verbose:
      sys.stderr.write(f"# processing skip list... ({skip_arg})\n")
    with open(skip_arg) as f:
      for line in f: skip.add(line.rstrip().split("\t")[0])
    if verbose:
      sys.stderr.write("# done: skipping computation for "+\
                       f"up to {len(skip)} units\n")
  elif verbose:
    sys.stderr.write("# no skip list, all input units will be processed")
  return skip

def get_mod_function(fn, fun, verbose):
  if fn:
    pmod = mod.py_or_nim(fn, verbose)
    return getattr(pmod, fun)
  else:
    return None

def open_files(args_out, args_log):
  outfile = open(args_out, "a") if args_out else sys.stdout
  logfile = open(args_log, "a") if args_log else sys.stderr
  return outfile, logfile

def close_files(outfile, logfile):
  if outfile != sys.stdout: outfile.close()
  if logfile != sys.stderr: logfile.close()

def compute_ids(unit_name, is_filename, idproc):
  identifier = idproc(unit_name) if idproc else unit_name
  return (unit_name if is_filename else identifier), identifier

def main(args):
  skip = compute_skip_set(args["--skip"], args["--verbose"])
  plugin = mod.py_or_nim(args["<plugin>"], args["--verbose"])
  idproc = get_mod_function(args["--idsproc"], "compute_id",
                            args["--verbose"])
  report = reports.Report.from_args(plugin, args)
  outfile, logfile = open_files(args["--out"], args["--log"])
  n_processed = 0
  params = args["--params"]
  state = plugin.initialize(**params.get("state", {})) \
      if hasattr(plugin, "initialize") else None
  for unit_name in tqdm.tqdm(input_units(args)):
    input_id, output_id = compute_ids(unit_name, args["<globpattern>"], idproc)
    if output_id in skip:
      skip.remove(output_id)
    else:
      compute(plugin, state, input_id, output_id, params,
              outfile, logfile, report, n_processed, args["--verbose"])
      n_processed += 1
  report.finalize(n_processed)
  if hasattr(plugin, "finalize"):
    plugin.finalize(state)
  close_files(outfile, logfile)

def validated(args):
  return scripts.validate(args, reports.args_schema,
      {"<globpattern>": Or(None, str), "<idsfile>": Or(None, os.path.exists),
       "--idsproc": Or(None, os.path.exists), "<col>": valid.optcolnum,
       "<plugin>": os.path.exists, "--out": Or(None, str),
       "--log": Or(None, str), "--skip": Or(None, os.path.exists)})

if "snakemake" in globals():
  args = snake.args(snakemake, reports.snake_args,
           input=["<plugin>", "--idsproc"],
           log=["--out", "--log"],
           params=["<globpattern>", "<idsfile>", "<col>", "--verbose",
                   "--skip"])
  if args["--skip"] is None:
    if args["--out"] and os.path.exists(args["--out"]):
      args["--skip"] = args["--out"]
    else:
      args["--skip"] = None
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(report_opts=reports.args_doc,
    common=scripts.args_doc), version="0.1")
  main(validated(args))
