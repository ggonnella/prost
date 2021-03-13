#!/usr/bin/env python3
"""
Compute GC percent from Fasta.
"""

import sh

statname = "gc_percent"
source = "computed from genomic.fas.gz; obtained by FTP from NCBI Refseq"

def value(fn):
  genomelen = int(str(sh.wc(sh.tr(sh.grep(sh.zcat(fn, _piped=True), ">", v=True,
      _piped=True), d="[:space:]", _piped=True), c=True)).rstrip())
  gclen = int(str(sh.wc(sh.tr(sh.sed(sh.grep(sh.zcat(fn, _piped=True), ">",
      v=True, _piped=True), r's/[^gcGC]//g', _piped=True), d="[:space:]",
      _piped=True), c=True)).rstrip())
  return str(gclen/genomelen)

