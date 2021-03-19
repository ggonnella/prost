"""
---
id: fas_stats_posix
version: 0.1.0
input: genome sequence; Fasta; optionally Gzip-compressed
output: genome_size, GC_content
method: count bases
implementation: pipe of posix tools using sh library
req_software: grep, tr, wc, zcat
parameters: gzipped (bool)
"""
