import sh
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

def nim(filename, verbose=False, import_constants=True):
    modulename = Path(filename).stem
    parent = Path(filename).parent
    m = nimporter.Nimporter.import_nim_module(modulename, [parent])
    if import_constants:
      for k in list(m.__dict__.keys()):
        if k.startswith("py_const_"):
          c = k[len("py_const_"):]
          setattr(m, c, getattr(m, k)())
          delattr(m, k)
    if verbose:
      sys.stderr.write(
          f"# nim module {modulename} imported from file {filename}\n")
    return m

def py_or_nim(filename, verbose=False):
  try:
    m = importer(filename, verbose)
    setattr(m, "__lang__", "nim")
  except:
    m = nim(filename, verbose)
    setattr(m, "__lang__", "python")
  return m

toml_content = """
[package]
authors = ["*"]
name = "{modulename}"
version = "1.0.0"
edition = "2018"

[dependencies]
pyo3 = {{ version = "0.13.2", features = ["abi3-py36", "extension-module"] }}
{more_dependencies}

[lib]
name = "{modulename}"
crate-type = ["cdylib"]
path = "{modulename}.rs"
"""

def rust(filename, verbose=False, import_constants=True):
  try:
    more_dependencies=sh.grep('\/\/ \[dependencies\] .*',
                              filename, P=True)
    more_dependencies=str(more_dependencies)[18:]
  except sh.ErrorReturnCode_1:
    more_dependencies=""
  workingdir = Path(filename).parent
  modulename = Path(filename).stem
  with open(workingdir/"Cargo.toml", "w") as f:
    f.write(toml_content.format(modulename = modulename,
                                more_dependencies=more_dependencies))
  sh.maturin("build", release=True, m=workingdir/"Cargo.toml")
  libfilename = workingdir/f"{modulename}.so"
  sh.ln(f"target/release/lib{modulename}.so", libfilename, s=True, f=True)
  sh.rm(workingdir/"Cargo.toml")
  sh.rm(workingdir/"Cargo.lock")
  spec = importlib.util.spec_from_file_location(modulename, libfilename)
  m = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(m)
  if import_constants:
    for k in list(m.Constants.__dict__.keys()):
      if not k.startswith("_"):
        setattr(m, k, getattr(m.Constants, k))
    delattr(m, "Constants")
  if verbose:
    sys.stderr.write(
        f"# python module {modulename} imported from file {filename}\n")
  return m

def py_or_nim_or_rust(filename, verbose=False):
  if Path(filename).suffix == ".rs":
    m = rust(filename, verbose)
    setattr(m, "__lang__", "rust")
  else:
    m = py_or_nim(filename, verbose)
  return m
