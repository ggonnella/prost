#!/usr/bin/env python3
"""
Count different types of features from a feature table
"""

counters = {
    "protein_coding_genes":
        lambda r: r["feature"] == "CDS" and \
                  r["class"]   == "with_protein",
    "rRNA_16S":
        lambda r: r["feature"] == "rRNA" and \
                  r["name"].startswith("16S"),
    "rRNA_5S":
        lambda r: r["feature"] == "rRNA" and \
                  r["name"].startswith("5S"),
    "rRNA_23S":
        lambda r: r["feature"] == "rRNA" and \
                  r["name"].startswith("23S"),
    "tRNA":
        lambda r: r["feature"] == "tRNA",
    "ncRNA":
        lambda r: r["feature"] == "ncRNA",
    "tmRNA":
        lambda r: r["feature"] == "tmRNA",
    "repeat_region":
        lambda r: r["feature"] == "repeat_region",
    "pseudogenes":
        lambda r: r["feature"] == "CDS" and \
                  r["class"] == "without_protein",
    "gene":
        lambda r: r["feature"] == "gene",
  }

def check(r):
  if r["feature"] not in ["tRNA", "ncRNA", "tmRNA", "repeat_region",
                          "gene", "CDS", "rRNA"]:
    sys.stderr.write(f"ERROR: unexpected feature type found: {r['feature']}\n")
  elif r["feature"] == "CDS" and \
       r["class"] not in ["with_protein", "without_protein"]:
    sys.stderr.write(f"ERROR: unexpected CDS class found: {r['class']}\n")
  elif r["feature"] == "rRNA" and \
       not r["name"].startswith("5S") and \
       not r["name"].startswith("16S") and \
       not r["name"].startswith("23S"):
    sys.stderr.write(f"ERROR: unexpected rRNA name found: {r['name']}\n")
