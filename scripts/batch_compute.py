#!/usr/bin/env python3
"""
Perform computations on multiple files, using the compute function of the
specified Python/Nim/Rust plugin module (see plugins/README.md for the plugins
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
  --idsproc FNAME   Python/Nim/Rust module, providing compute_id(str)->str;
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
  --serial          run the computation serially (default: use multiprocessing)
{report_opts}
{common}
"""

from docopt import docopt
from schema import Or
import os
import sys
from pathlib import Path
from glob import glob
from lib import snake, valid, reports, scripts, formatting, plugins
import tqdm
from concurrent.futures import as_completed, ProcessPoolExecutor
from functools import partial
import multiplug

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
    sys.stderr.write("# no skip list, all input units will be processed\n")
  return skip

def get_mod_function(filename, fun, verbose):
  if filename:
    pmod = multiplug.importer(filename, verbose=verbose,
                              **plugins.IDPROC_PLUGIN_INTERFACE)
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

def compute_all_ids(args, idproc, skip, verbose):
  if verbose:
    sys.stderr.write("# compute input and output IDs...")
  result = []
  for unit_name in input_units(args):
    input_id, output_id = compute_ids(unit_name, args["<globpattern>"], idproc)
    if output_id in skip:
      skip.remove(output_id)
    else:
      result.append((input_id, output_id))
  if verbose:
    sys.stderr.write("done\n")
  return result

def init_state_and_get_params(params, plugin):
  if plugin.initialize is not None:
    params["state"] = plugin.initialize(**params.get("state", {}))
  elif "state" not in params:
    params["state"] = None
  return params

def on_failure(outfile, logfile, report, output_id, exc):
  outfile.flush()
  logfile.flush()
  report.error(exc, output_id)

def on_success(outfile, logfile, report, output_id, results, logs):
  results = "\t".join([str(r) for r in results])
  outfile.write(f"{output_id}\t{results}\n")
  if logfile and logs:
    for msg in logs:
      logfile.write(f"{output_id}\t{msg}\n")
  report.step()

plugin = None

def unit_processor(input_id, params):
  global plugin
  return plugin.compute(input_id, **params)

def run_in_parallel(all_ids, params, outfile, logfile, report, desc, verbose):
  if verbose:
    sys.stderr.write("# Computation will be in parallel (multiprocess)\n")
  with ProcessPoolExecutor() as executor:
    futures_map = {executor.submit(unit_processor, unit_ids[0], params):
                   unit_ids[1] for unit_ids in all_ids}
    for future in tqdm.tqdm(as_completed(futures_map), total=len(futures_map),
                            desc=desc):
      output_id = futures_map[future]
      try:
        results, logs = future.result()
      except Exception as exc:
        on_failure(outfile, logfile, report, output_id, exc)
        raise(exc)
      else:
        on_success(outfile, logfile, report, output_id, results, logs)

def run_serially(all_ids, params, outfile, logfile, report, desc, verbose):
  if verbose:
    sys.stderr.write("# Computation will be serial\n")
  global plugin
  for unit_ids in tqdm.tqdm(all_ids, desc=desc):
    output_id = unit_ids[1]
    try:
      results, logs = plugin.compute(unit_ids[0], **params)
    except Exception as exc:
      on_failure(outfile, logfile, report, output_id, exc)
      raise(exc)
    else:
      on_success(outfile, logfile, report, output_id, results, logs)

def main(args):
  global plugin
  skip = compute_skip_set(args["--skip"], args["--verbose"])
  plugin = multiplug.importer(args["<plugin>"], verbose=args["--verbose"],
                              **plugins.COMPUTE_PLUGIN_INTERFACE)
  desc = formatting.shorten(Path(args["<plugin>"]).stem, 15)
  idproc = get_mod_function(args["--idsproc"], "compute_id",
                            args["--verbose"])
  report = reports.Report.from_args(plugin, args)
  outfile, logfile = open_files(args["--out"], args["--log"])
  params = init_state_and_get_params(args["--params"], plugin)
  all_ids = compute_all_ids(args, idproc, skip, args["--verbose"])
  run = run_serially if args["--serial"] else run_in_parallel
  run(all_ids, params, outfile, logfile, report, desc, args["--verbose"])
  report.finalize()
  if plugin.finalize is not None:
    plugin.finalize(params["state"])
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
                   "--skip", "--serial"])
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
