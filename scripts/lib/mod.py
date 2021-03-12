import sys
import importlib
from pathlib import Path

def importer(filename, verbose=False):
  modulename = Path(filename).stem
  spec = importlib.util.spec_from_file_location(modulename, filename)
  m = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(m)
  if verbose:
    sys.stderr.write(f"# module {modulename} imported from file {filename}\n")
  return m
