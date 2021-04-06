#!/usr/bin/env python3
"""
Update the table of Genbank assembly accessions to Refseq assembly accessions
and Genbank sequence accessions.

Usage:
  gb_seq_accessions_for_assemblies.py [options] <summary>

Arguments:
  summary: Refseq assembly_summary.txt file of NCBI

Options:
  --complete, -c            only "Complete Genome" and marked as "latest"
  --update, -u F            file with previous results, to update
{common}
"""

from docopt import docopt
import sh
import os
import sys
import tqdm
from lib import snake, scripts
from schema import Or

def is_complete(elems):
  """
  Filter for the lines of a assembly_summary table.
  True if the latest available assembly of a complete genome.
  """
  return elems[10] == "latest" and elems[11] == "Complete Genome"

def gb_ftppath_and_filename(gbasm, rspath, sfx = "_genomic.fna.gz"):
  """
  Given a FTP path for the Refseq assembly and the equivalent
  Genbank accession, compute the Genbank FTP path and filename.
  """
  # e.g. ['ftp:', '', 'ftp.ncbi.nlm.nih.gov', 'genomes', 'all',
  #       'GCF', '001', '458', '655', 'GCF_001458655.1_X_Y_Z']
  path_parts = rspath.split("/")
  # extract the accession parts which make up directory names in the final path
  # e.g. "GCA_002287175.1" => "GCA", "002", "287", "175"
  subdirnames = [gbasm[0:3], gbasm[4:7], gbasm[7:10], gbasm[10:13]]
  path_parts[5:9] = subdirnames
  # extract the assembly name from the last part of the path (e.g. X_Y_Z)
  # and replace the assembly accession with the Genbank one
  asmname = path_parts[9].split("_", 2)[-1]
  asmdir = gbasm + "_" + asmname
  path_parts[9] = asmdir
  gbftppath = "/".join(path_parts)
  gbfilename = asmdir + sfx
  return gbftppath, gbfilename

def download_and_extract_fastaids(remotefn, localfn, maxtries=3):
  i = 0
  while i < maxtries:
    try:
      sh.wget(remotefn, O=localfn)
      # zcat localfn | grep -P '^>' | cut -f1 -d' ' | cut -c2-
      seqaccs=sh.cut(sh.cut(sh.grep(sh.zcat(localfn,_piped=True),
        "^>",P=True,_piped=True), f="1",d=" ",_piped=True),
        c="2-").split("\n")[:-1]
      os.remove(localfn)
      return seqaccs
    except sh.ErrorReturnCode:
      i+=1
      continue
  if os.path.exists(localfn):
    os.remove(localfn)
  return None

def compute_known(filename):
  result = set()
  if filename and os.path.exists(filename):
    with open(filename) as f:
      for line in f:
        elems = line.rstrip().split("\t")
        result.add(elems[0])
  return result

def compute_targets(summaryfile, known_results, completeonly):
  result = []
  with open(summaryfile) as f:
    for line in f:
      if line[0] == "#": continue
      elems = line.rstrip().split("\t")
      rsasm = elems[0]
      gbasm = elems[17]
      rspath = elems[19]
      if rspath == "na": continue
      if completeonly and not is_complete(elems): continue
      if gbasm in known_results: continue
      result.append((gbasm, rspath, rsasm))
  return result

def compute_accessions(targets, outfile):
  noferrors = 0
  for gbasm, rspath, rsasm in tqdm.tqdm(targets):
    ftppath, fn = gb_ftppath_and_filename(gbasm, rspath)
    sys.stderr.write(f"Downloading {fn}...\n")
    seqaccs = download_and_extract_fastaids(os.path.join(ftppath, fn), fn)
    if seqaccs is None:
      sys.stderr.write(f"ERROR: it was not possible to download {fn}\n")
      noferrors += 1
    outfile.write("\t".join([gbasm, rsasm, ",".join(seqaccs)])+"\n")
    outfile.flush()
  sys.stderr.write(f"Number of errors: {noferrors}\n")

def open_outfile(args):
  if args["--update"]:
    return open(args["--update"], "a")
  else:
    return sys.stdout

def close_outfile(outfile, args):
  if args["--update"]:
    outfile.close()

def main(args):
  known_results = compute_known(args["--update"])
  targets = compute_targets(args["<summary>"], known_results,
                            args["--complete"])
  outfile = open_outfile(args)
  compute_accessions(targets, outfile)
  close_outfile(outfile, args)

def validated(args):
  return scripts.validate(args, {"--update": Or(None, str),
    "<summary>": open, "--complete": Or(None, bool)})

if "snakemake" in globals():
  args = snake.args(snakemake, params=[("--update", "prev")],
                    input=["<summary>"])
  args["--complete"]=True
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(common=scripts.args_doc), version=0.1)
  main(validated(args))

