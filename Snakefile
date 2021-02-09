#!/usr/bin/env python3
"""

Task 2: update_accession_tables
  2a: refseq
  - extract the accession numbers -- for not yet done files; and update the tables
  2b: genbank
  - temporarily download also the genbank for not yet done files, to update the accession tables

Task 3: update_bacdive
  update the bacdive database

Task 4: associate_complete_genomes_to_bacdive
  4a:
  - match the bacdive records using the accession tables
  4b:
  - NEW: match the genomes using the IDs

Task 5: analyse_bacdive_details
  - structure;
  - NEW: single values

Task 6: select_genomes_for_bacdive_trait
  - collect the genomes corresponding to a value / all values of a trait
"""

from glob import glob
import os
import shutil
import re

DBs = ["refseq", "genbank"]
Domains = ["bacteria", "archaea"]

genomesdir=os.environ["HOME"]+"/data/prok_complete_genomes/"
asmsummarydir=genomesdir+"asmsummary/"

#onstart:
#  for db in DBs:
#    for domain in Domains:
#      output=asmsummarydir+f"assembly_summary_{db}_{domain}.txt"
#      remote=f"ftp://ftp.ncbi.nlm.nih.gov/genomes/{db}/{domain}/assembly_summary.txt"
#      shell(f"curl -o {output} -z {output} {remote}")

rule all:
  input:
    expand(asmsummarydir+"plan/complete_genomes_refseq_{domain}.{unit}.done",
        domain=Domains,unit=["new", "obsolete"])

rule list_complete_genomes_for_domain:
  input:
    asmsummarydir+"assembly_summary_{db}_{domain}.txt"
  output:
    asmsummarydir+"complete_genomes_{db}_{domain}.txt"
  run:
    o_f = open(output[0], "w")
    with open(input[0]) as f:
      for line in f:
        if line[0] != "#":
          elems = line.rstrip().split("\t")
          if elems[10] == "latest" and elems[11] == "Complete Genome":
            # asmacc ftppath
            o_f.write("\t".join([elems[0], elems[19]])+"\n")
    o_f.close()

rule create_updating_plan:
  input:
    asmsummarydir+"complete_genomes_refseq_{domain}.txt"
  output:
    new=asmsummarydir+"plan/complete_genomes_refseq_{domain}.new",
    upt=asmsummarydir+"plan/complete_genomes_refseq_{domain}.uptodate",
    obs=asmsummarydir+"plan/complete_genomes_refseq_{domain}.obsolete"
  params:
    chunkspfx=genomesdir+"seq/{domain}/chunk."
  run:
    filenames=set(os.path.basename(fn) for fn in \
        glob(params.chunkspfx+"*/*_genomic.fna.gz"))
    out_new = open(output.new, "w")
    out_obs = open(output.obs, "w")
    out_upt = open(output.upt, "w")
    with open(input[0]) as f:
      for line in f:
        elems = line.rstrip().split("\t")
        fn = elems[1].split("/")[-1] + "_genomic.fna.gz"
        if fn in filenames:
          out_upt.write(fn+"\n")
          filenames.remove(fn)
        else:
          out_new.write(fn+"\t"+elems[1]+"\n")
      for fn in filenames:
        out_obs.write(fn+"\n")
    out_new.close()
    out_obs.close()
    out_upt.close()

rule move_obsolete_genomes:
  input:
    asmsummarydir+"plan/complete_genomes_refseq_{domain}.obsolete"
  output:
    touch(asmsummarydir+"plan/complete_genomes_refseq_{domain}.obsolete.done")
  params:
    chunkspfx=genomesdir+"seq/{domain}/chunk."
  run:
    with open(input[0]) as f:
      for line in f:
        fn = line.rstrip()
        srcfile = glob(params.chunkspfx+"*/"+fn)[0]
        obsdir = os.path.dirname(srcfile)+"/obsolete"
        if not os.path.exists(obsdir):
          os.mkdir(obsdir)
        shutil.move(srcfile, obsdir)

rule download_new_genomes:
  input:
    asmsummarydir+"plan/complete_genomes_refseq_{domain}.new"
  output:
    touch(asmsummarydir+"plan/complete_genomes_refseq_{domain}.new.done")
  params:
    chunkspfx=genomesdir+"seq/{domain}/chunk.",
    chunkmaxsize=500
  run:
    if os.stat(input[0]).st_size > 0:
      chunknums=[int(os.path.splitext(x)[1][1:]) for x in glob(params.chunkspfx+"*")]
      if len(chunknums) == 0:
        nextchunk= 0
      else:
        nextchunk=max(chunknums)+1
      nextchunkdir=params.chunkspfx+str(nextchunk)
      if not os.path.exists(nextchunkdir):
        os.mkdir(nextchunkdir)
      genomenum=0
      with open(input[0]) as f:
        for line in f:
          elems=line.rstrip().split("\t")
          if genomenum == params.chunkmaxsize:
            genomenum=0
            nextchunk+=1
            nextchunkdir=params.chunkspfx+str(nextchunk)
            os.mkdir(nextchunkdir)
          genomenum+=1
          shell("wget "+elems[1]+"/"+elems[0]+" -O "+nextchunkdir+"/"+elems[0])

rule table_asm_accessions:
  input:
    asmsummarydir+"assembly_summary_refseq_{domain}.txt"
  output:
    asmsummarydir+"accessions/rsasm.{domain}.tsv",
    asmsummarydir+"accessions/gbasm_rsasm.{domain}.tsv"
  run:
    r_tsv = open(output[0], "w")
    g_tsv = open(output[1], "w")
    with open(input[0]) as f:
      for line in f:
        if line[0] != "#":
          elems = line.rstrip().split("\t")
          if elems[10] == "latest" and elems[11] == "Complete Genome":
            r_tsv.write(elems[0]+"\n")
            # gbasmacc asmacc
            g_tsv.write(elems[17]+"\t"+elems[0]+"\n")
    r_tsv.close()
    g_tsv.close()

rule genbank_download_paths:
  input:
    asmsummarydir+"accessions/gbasm_rsasm.{domain}.tsv",
    asmsummarydir+"assembly_summary_genbank_{domain}.txt"
  output:
    asmsummarydir+"accessions/plan/genbank_{domain}.remote_path.tsv"
  run:
    gbasm = dict()
    outfile = open(output[0], "w")
    with open(input[0]) as f:
      for line in f:
        elems = line.rstrip().split("\t")
        gbasm[elems[0]] = elems[1]
    with open(input[1]) as f:
      for line in f:
        if line[0] != "#":
          elems = line.rstrip().split("\t")
          if elems[0] in gbasm:
            rsasm = gbasm[elems[0]]
            del gbasm[elems[0]]
            outfile.write("\t".join([rsasm, elems[0], elems[19]])+"\n")
    outfile.close()

rule cp_existing_genbank_seqacc:
  params:
    possible_input=asmsummarydir+"accessions/gbasm_rsasm_gbseq.{domain}.tsv"
  output:
    asmsummarydir+"accessions/plan/gbasm_rsasm_gbseq.{domain}.tsv"
  shell:
    """
    if [ -e {params.possible_input} ]; then
      cp {params.possible_input} {output}
    else
      touch {output}
    fi
    """

rule list_needed_genbank_seqacc:
  input:
    asmsummarydir+"accessions/plan/genbank_{domain}.remote_path.tsv",
    asmsummarydir+"accessions/plan/gbasm_rsasm_gbseq.{domain}.tsv"
  output:
    asmsummarydir+"accessions/plan/genbank_{domain}.need_seqacc"
  run:
    not_needed = set()
    with open(input[1]) as f:
      for line in f:
        elems = line.rstrip().split("\t")
        not_needed.add(elems[0])
    outfile = open(output[0], "w")
    with open(input[0]) as f:
      for line in f:
        elems = line.rstrip().split("\t")
        if elems[1] not in not_needed:
          outfile.write(line)
    outfile.close()

rule update_gbseq_list:
  input:
    asmsummarydir+"accessions/plan/genbank_{domain}.need_seqacc",
    asmsummarydir+"accessions/plan/gbasm_rsasm_gbseq.{domain}.tsv"
  output:
    asmsummarydir+"accessions/gbasm_rsasm_gbseq.{domain}.tsv"
  run:
    outfile = open(output[0], "w")
    with open(input[0]) as f:
      for line in f:
        elems = line.rstrip().split("\t")
        fpath = elems[2]
        fn = elems[2].split("/")[-1] + "_genomic.fna.gz"
        accfn = elems[1] + ".acc"
        shell(f"wget {fpath}/{fn} -O {fn}")
        shell(f"zcat {fn} | grep -P '^>' | cut -f1 -d' ' | cut -c2- > {accfn}")
        seqaccs = []
        with open(accfn) as f2:
          for l2 in f2:
            seqaccs.append(l2.rstrip())
        outfile.write("\t".join([elems[1], elems[0], ",".join(seqaccs)])+"\n")
        outfile.flush()
        os.remove(fn)
        os.remove(accfn)
    with open(input[1]) as f:
      for line in f:
        outfile.write(line)
    outfile.close()
