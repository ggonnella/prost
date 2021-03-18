#!/usr/bin/env python3
"""
Compute genome length from Fasta.

---
id:              fas_stats_nim
version:         0.1.0
input:           genome sequence; Fasta; optionally Gzip-compressed
output:          genome_size, GC_content
method:          count bases
implementation:  NIM module, lookup table
parameters:      gzipped (Boolean)
req_software:    nim >= 1.2.8; nimble:zip >= 0.3.1
req_hardware:    NULL
advice:          preferred to fasgz_stats_posix, if nim is available,
                 since it is faster
"""

import sh

statname = "genome_length"
source = "computed from genomic.fas.gz; obtained by FTP from NCBI Refseq"

def value(fn):
  return str(sh.wc(sh.tr(sh.grep(sh.zcat(fn, _piped=True), ">", v=True,
         _piped=True), d="[:space:]", _piped=True), c=True)).rstrip()

