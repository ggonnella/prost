#!/usr/bin/env python3
"""
Count features of different type and class
from a NCBI Refseq genome annotation GFF file
"""

from collections import defaultdict

def results(db):
  counters = defaultdict(int)
  genes = db.features_of_type("gene")
  for g in genes:
    children = list(db.children(g.id, level=1))
    assert(len(children) >= 1)
    child = children[0]
    assert(all([child.featuretype == c.featuretype for c in children]))
    if child.featuretype == "CDS":
      counters["protein_coding_genes"] += 1
    elif child.featuretype == "tRNA":
      counters["tRNAs"] += 1
    elif child.featuretype in ["ncRNA", "RNase_P_RNA", "SRP_RNA", "tmRNA",
                               "antisense_RNA"]:
      counters["ncRNAs"] += 1
    elif child.featuretype == "rRNA":
      product = child.attributes['product'][0]
      if product.startswith("16S"):
        counters["rRNA_16S"] += 1
      elif product.startswith("5S"):
        counters["rRNA_5S"] += 1
      elif product.startswith("23S"):
        counters["rRNA_23S"] += 1
      else:
        assert(False)
    else:
      raise ValueError(f"Feature Type under Gene: {child.featuretype} "+\
                       f"unexpected (ID={child.id})")
  counters["pseudogenes"] = db.count_features_of_type("pseudogene")
  return counters

