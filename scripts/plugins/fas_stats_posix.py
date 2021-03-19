"""
Compute sequence statistics using Posix command line tools
"""
import sh

def analyze(filename, **kwargs):
  """
  # Plugin documentation in YAML format:
  id: fas_stats_posix
  version: 0.1.0
  input: genome sequence; Fasta; optionally Gzip-compressed
  output: sequence_attribute.(genome_size,GC_content)
  method: count bases
  implementation: pipe of posix tools using sh library
  req_software: grep, tr, wc, zcat
  parameters: uncompressed (bool)
  """
  results = {}
  countchars = sh.wc.bake(c=True)
  gc_only = sh.tr.bake("-dc", "[GCgc]", _piped=True)
  non_space = sh.tr.bake(d="[:space:]", _piped=True)
  vgrep = sh.grep.bake(_piped=True, v=True)
  uncompressed = sh.zcat.bake(_piped=True)
  if kwargs.get("uncompressed", False):
    genomelen = int(countchars(non_space(vgrep("^>", filename))))
    gclen = int(countchars(gc_only(vgrep("^>", filename))))
  else:
    genomelen = int(countchars(non_space(vgrep(uncompressed(filename), "^>"))))
    gclen = int(countchars(gc_only(vgrep(uncompressed(filename), "^>"))))
  results["genome_size"] = genomelen
  results["gc_content"] = gclen/genomelen
  return {"sequence_attribute": results}, {}
