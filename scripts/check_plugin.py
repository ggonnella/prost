#!/usr/bin/env python3
"""
Check that a computation plugin respects the specification
given in plugins/README.md.

The following is checked:
- the signature of the compute function and, optionally of
  the initialize, finalize (for Python plugins only)
- the existance, type and value of the plugin constants
- all the attributes declared in the OUTPUT constant must be
  contained in the given attribute definitions file

Usage:
  check_plugin.py <plugin> <definitions> [options]

Arguments:
  plugin:       Python or Nim batch_compute.py plugin module
  definitions:  YAML file containing the attribute definitions

Options:
{common}
"""

from schema import And, Use
import yaml
import sys
import os
import inspect
from lib import scripts
import multiplug
import snacli

def main(args):
  plugin = multiplug.importer(args["<plugin>"], verbose=args["--verbose"])
  exit_code = 0
  if not hasattr(plugin, "compute"):
    sys.stderr.write("# [Error] plugin does not provide a compute function\n")
    exit_code = 1
  else:
    if args["--verbose"]:
      sys.stderr.write("# [OK] plugin provides a compute function\n")
    if plugin.__lang__ == "python":
      compute_spec = inspect.getfullargspec(plugin.compute)
      if compute_spec.varargs is not None:
        sys.stderr.write("# [Error] "+\
            "plugin.compute() accepts variable positional arguments\n")
        exit_code = 1
      if not hasattr(plugin, "initialize"):
        if len(compute_spec.args) != 1:
          sys.stderr.write("# [Error] "+\
              "plugin.compute() does not accept a single positional argument\n")
          exit_code = 1
      else:
        if (len(compute_spec.args) != 2) or (compute_spec.args[1] != "state"):
          sys.stderr.write("# [Error] "+\
              "plugin.compute() in plugin with initialize function "+\
              "does not accept a state keyword argument\n")
          exit_code = 1
      if compute_spec.varkw is None:
        sys.stderr.write("# [Error] "+\
            "plugin.compute() does not accept variable keywords arguments\n")
        exit_code = 1
      if hasattr(plugin, "initialize"):
        init_spec = inspect.getfullargspec(plugin.initialize)
        if (init_spec.varargs is not None) or (len(init_spec.args) > 0):
          sys.stderr.write("# [Error] "+\
              "plugin.initialize() accepts positional arguments\n")
          exit_code = 1
        if init_spec.varkw is None:
          sys.stderr.write("# [Error] "+\
              "plugin.initialize() does not accept variable keywords arguments\n")
          exit_code = 1
      if hasattr(plugin, "finalize"):
        f_spec = inspect.getfullargspec(plugin.finalize)
        if f_spec.varargs is not None:
          sys.stderr.write("# [Error] "+\
              "plugin.finalize() accepts variable positional arguments\n")
          exit_code = 1
        if len(f_spec.args) != 1:
          sys.stderr.write("# [Error] "+\
              "plugin.finalize() does not accept a single positional argument\n")
          exit_code = 1
        if f_spec.varkw is not None:
          sys.stderr.write("# [Error] "+\
              "plugin.initialize() accepts variable keywords arguments\n")
          exit_code = 1
    elif args["--verbose"]:
      sys.stderr.write("# [Info] signature of plugin functions not analyzed, "+\
          " since it is a {} plugin\n", plugin.__lang__)
  if hasattr(plugin, "finalize"):
    if not hasattr(plugin, "initialize"):
      sys.stderr.write("# [Error] "+\
            "plugin.finalize() defined in plugin without plugin.initialize()\n")
      exit_code = 1
  for const in ["ID", "VERSION", "INPUT", "OUTPUT"]:
    if not hasattr(plugin, const):
      sys.stderr.write("# [Error] "+\
          f"plugin does not define a constant named '{const}'\n")
      exit_code = 1
    elif args["--verbose"]:
      sys.stderr.write(f"[OK] plugin defines a constant named '{const}'\n")
  for const, maxlen in {"ID": 256, "VERSION": 64, "INPUT": 512,
                        "METHOD": 4096, "IMPLEMENTATION": 4096,
                        "REQ_SOFTWARE": 4096, "REQ_HARDWARE": 4096,
                        "ADVICE": 4096}.items():
    if hasattr(plugin, const):
      v = getattr(plugin, const)
      if v is not None:
        if not isinstance(v, str):
          sys.stderr.write(f"# [Error] plugin.{const} must be a string\n")
          exit_code = 1
        elif len(v) > maxlen:
          sys.stderr.write(f"# [Error] plugin.{const} length must <= {maxlen}\n")
          exit_code = 1
        elif args["--verbose"]:
          sys.stderr.write(f"# [OK] plugin.{const} type and format is valid\n")
  if hasattr(plugin, "OUTPUT"):
    if not isinstance(plugin.OUTPUT, list):
      sys.stderr.write("# [ERROR] plugin.OUTPUT is not a list\n")
      exit_code = 1
    notstr = [key for key in plugin.OUTPUT if not isinstance(key, str)]
    if notstr:
      sys.stderr.write("# [Error] "+\
          f"plugin.OUTPUT non-string elements found: {notstr}\n")
      exit_code = 1
    elif args["--verbose"]:
      sys.stderr.write("# [OK] all plugin.OUTPUT elements are strings\n")
    notin = [key for key in plugin.OUTPUT \
              if key not in notstr and key not in args["<definitions>"]]
    if notin:
      sys.stderr.write("# [Error] plugin.OUTPUT elements "+\
          f"not found in the attributes definitions file: {notin}\n")
      exit_code = 1
    elif args["--verbose"]:
      sys.stderr.write("# [OK] all plugin.OUTPUT elements are "+\
            "present in the attributes definitions file\n")
  if hasattr(plugin, "PARAMETERS") and plugin.PARAMETERS is not None:
    if not isinstance(plugin.PARAMETERS, list):
      sys.stderr.write("# [ERROR] plugin.PARAMETERS is not a list\n")
      exit_code = 1
    for element in plugin.PARAMETERS:
      if not isinstance(element, tuple) and not isinstance(element, list):
        sys.stderr.write("# [ERROR] "+\
            f"plugin.PARAMETER element not tuple: {element}\n")
        exit_code = 1
      elif len(element) != 4:
        sys.stderr.write("# [ERROR] "+\
            f"plugin.PARAMETER element tuple length is not 4: {element}\n")
        exit_code = 1
      else:
        for s in element:
          if not isinstance(s, str):
            sys.stderr.write("# [ERROR] "+\
              "plugin.PARAMETER tuple element is not a "+\
              f"string: {s} of {element}\n")
            exit_code = 1

  return exit_code

def validated(args):
  return scripts.validate(args, {"<plugin>": os.path.exists,
    "<definitions>": And(str, Use(open), Use(yaml.safe_load))})

with snacli.args(input=["<plugin>", "<definitions>"],
                 docvars={"common": scripts.args_doc},
                 version="0.1") as args:
  if args: main(validated(args))
