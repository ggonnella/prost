import sys
import importlib
from pathlib import Path
import nimporter

def importer(filename, verbose=False):
  modulename = Path(filename).stem
  spec = importlib.util.spec_from_file_location(modulename, filename)
  m = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(m)
  if verbose:
    sys.stderr.write(
        f"# python module {modulename} imported from file {filename}\n")
  return m

def nim(filename, verbose=False):
    modulename = Path(filename).stem
    parent = Path(filename).parent
    m = nimporter.Nimporter.import_nim_module(modulename, [parent])
    if verbose:
      sys.stderr.write(
          f"# nim module {modulename} imported from file {filename}\n")
    return m

def py_or_nim(filename, verbose=False):
  try:
    m = importer(filename, verbose)
    setattr(m, "__nim__", False)
  except:
    m = nim(filename, verbose)
    setattr(m, "__nim__", True)
  return m
