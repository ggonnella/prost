#!/usr/bin/env python3
"""
Compute GC percent from Fasta.
"""

import sh
import os
import nimporter, computestats.fas_stats as fas_stats

statname = "gc_percent"
source = "computed from genomic.fas.gz; obtained by FTP from NCBI Refseq"

def value_wo_zip_library(fn):
  sh.gunzip(fn, k=True)
  unzipped_fn = fn[:-3]
  result = fas_stats.gcpercent(unzipped_fn, False)
  os.remove(unzipped_fn)
  return str(result)

def value(fn):
  result = fas_stats.gcpercent(fn, True)
  return str(result)

