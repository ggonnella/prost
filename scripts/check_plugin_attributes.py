#!/usr/bin/env python3
"""
Check that all the attributes declared in the OUTPUT variable
of a plugin are already in an attribute definitions file.

Usage:
  check_plugin_attributes.py <plugin> <definitions> [options]

Arguments:
  plugin:       Python or Nim batch_compute.py plugin module
  definitions:  YAML file containing the attribute definitions

Options:
{common}
"""

from docopt import docopt
from schema import Or, And, Use
import yaml
import sys
import os
from lib import snake, mod, valid, scripts

def main(args):
  plugin = mod.py_or_nim(args["<plugin>"], args["--verbose"])
  notin = [key for key in plugin.OUTPUT if key not in args["<definitions>"]]
  if notin:
    if args["--verbose"]:
      sys.stderr.write("# Keys not found in the attributes definitions:\n")
    for key in notin:
      sys.stderr.write(f"{key}\n")
    return 1
  else:
    if args["--verbose"]:
      sys.stderr.write("# All keys found in the attributes definitions\n")
    return 0

def validated(args):
  return scripts.validate(args, {"<plugin>": os.path.exists,
    "<definitions>": And(str, Use(open), Use(yaml.safe_load))})

if "snakemake" in globals():
  args = snake.args(snakemake, input=["<plugin>", "<definitions>"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(common=scripts.args_doc), version="0.1")
  main(validated(args))
