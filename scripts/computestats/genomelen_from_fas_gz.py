#!/usr/bin/env python3
"""
Compute genome length from Fasta.
"""

import sh

statname = "genome_length"
source = "computed from genomic.fas.gz; obtained by FTP from NCBI Refseq"

def value(fn):
  return str(sh.wc(sh.tr(sh.grep(sh.zcat(fn, _piped=True), ">", v=True,
         _piped=True), d="[:space:]", _piped=True), c=True)).rstrip()

