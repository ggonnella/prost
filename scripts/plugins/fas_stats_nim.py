"""
---
id: fas_stats_nim
version: 0.1.0
input: genome sequence; Fasta; optionally Gzip-compressed
output: sequence_attribute.(genome_size,GC_content)
method: count bases
implementation: NIM module, lookup table
parameters: gzipped (bool)
req_software: nim >= 1.2.8; nimble:zip >= 0.3.1
advice: >
  preferred to fasgz_stats_posix, if nim is available, since it is faster
"""
