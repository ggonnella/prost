"""
Helper functions for scripts
"""

def setargs(args, src, *names):
  """
  For each name in <names>, set args[name] to
  src.get(key), where key is a transformation of the name
  as decribed below.

  This is to facilitate scripts to be called directly
  with option parsing by docopt or via snakemake "script".

  The <names> must be strings such those returned by docopt,
  i.e. starting with "--" or enclosed in "<>".
  The key is derived from the name by removing these formatting chars
  and replacing internal "-" with underscores.
  If a key is not available, args[name] is set to None.

  E.g. setargs(args, snakemake.input, "<file>")
       => args["<file>"] = snakemake.input.get("file")
       setargs(args, snakemake.params, "--long-opt")
       => args["--long-opt"] = snakemake.params.get("long_opt")
  """
  for name in names:
    if name.startswith("--"):
      key = name[2:]
    else:
      assert(name[0]=="<")
      assert(name[-1]==">")
      key = name[1:-1]
    key = key.replace("-","_")
    args[name] = src.get(key)

def args(snakemake, **kwargs):
  args = {}
  for k, v in kwargs.items():
    setargs(args, getattr(snakemake, k), *v)
  return args

