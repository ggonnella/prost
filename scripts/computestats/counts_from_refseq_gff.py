#!/usr/bin/env python3
"""
Count features of different type and class
from a NCBI Refseq genome annotation GFF file.
"""

from collections import defaultdict

def process_gene_feature(db, gene, counters, info):
  children = list(db.children(gene.id, level=1))
  if len(children) == 0:
    info["gene_with_no_children"].append(str(gene))
    return
  child = children[0]
  if not all([child.featuretype == c.featuretype for c in children]):
    info["gene_with_heterogeneous_children_type"].append(str(gene))
    return
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
      info["rRNA_unknown"].append(str(child))
      assert(False)
  elif child.featuretype == "hammerhead_ribozyme":
    counters["hammerhead_ribozymes"] += 1
  else:
    info["unexpected_gene_child_feature_type"].append(str(child))

def analyze(db):
  counters = defaultdict(int)
  info = defaultdict(list)
  for ftr in db.all_features():
    if len(list(db.parents(ftr.id))) == 0:
      if ftr.featuretype == "gene":
        process_gene_feature(db, ftr, counters, info)
      elif ftr.featuretype == "pseudogene":
        counters["pseudogenes"] += 1
      elif ftr.featuretype == "region":
        counters["n_molecules"] += 1
      elif ftr.featuretype in ["sequence_feature"]:
        info["generic_sequence_feature"].append(str(ftr))
      elif ftr.featuretype == "riboswitch":
        counters["riboswitches"] += 1
      elif ftr.featuretype == "direct_repeat":
        if ftr.attributes["rpt_family"][0] == "CRISPR":
          counters["CRISPRs"] += 1
        else:
          info["direct_repeat_not_CRISPR"].append(str(ftr))
      else:
        info["top_level_feature_other_type"].append(str(ftr))
  return counters, info

