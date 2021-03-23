#!/usr/bin/env python3
"""
Count features of different type and class
from a NCBI Refseq genome annotation GFF file.
"""

ID =      "refseq_gff_count"
VERSION = "0.1.0"
INPUT =   "genome annotation; from NCBI Refseq; gff, optionally gzipped"
OUTPUT =  [
  "n_protein_coding_gene", "n_tRNA", "n_ncRNA", "n_rRNA_16S", "n_rRNA_5S",
  "n_rRNA_23S", "n_hammerhead_ribozyme", "n_pseudogene", "n_molecule",
  "n_riboswitch", "n_CRISPR",
]
METHOD = """
Counts features of a predefined set of feature_type
- feature_type is defined as in GFF (i.e. term from the SOFA
  subset of the SO ontology)
- to avoid counting portions of the same feature multiple
  times, children of the same type of the same gene
  are counted only once
- for some feature_type (e.g. rRNA), a more specific
  feature_type (e.g. rRNA_16S) is assigned using information
  from the attributes (e.g. product_name)
"""
IMPLEMENTATION = "Based on the py library Gffutils, using an in-memory database"
REQ_SOFTWARE = "pip library gffutils"

from collections import defaultdict
import gffutils
import sys

def _dict2list(d, keynames):
  assert(not(set(d.keys())-set(keynames)))
  [d.get(k, 0) for k in keynames]

def _process_gene_feature(db, gene, counters, logs):
  children = list(db.children(gene.id, level=1))
  if len(children) == 0:
    logs.append(f"gene_with_no_children\t{gene}")
    return
  child = children[0]
  if not all([child.featuretype == c.featuretype for c in children]):
    logs.append(f"gene_with_heterogeneous_children_type\t{gene}")
    return
  if child.featuretype == "CDS":
    counters["protein_coding_gene"] += 1
  elif child.featuretype == "tRNA":
    counters["tRNA"] += 1
  elif child.featuretype in ["ncRNA", "RNase_P_RNA", "SRP_RNA", "tmRNA",
                             "antisense_RNA"]:
    counters["ncRNA"] += 1
  elif child.featuretype == "rRNA":
    product = child.attributes['product'][0]
    if product.startswith("16S"):
      counters["rRNA_16S"] += 1
    elif product.startswith("5S"):
      counters["rRNA_5S"] += 1
    elif product.startswith("23S"):
      counters["rRNA_23S"] += 1
    else:
      logs.append(f"rRNA_unknown\t{child}")
      assert(False)
  elif child.featuretype == "hammerhead_ribozyme":
    counters["hammerhead_ribozyme"] += 1
  else:
    logs.append(f"unexpected_gene_child_feature_type\t{child}")

def compute(filename, **kwargs):
  db = gffutils.create_db(filename, ":memory:",
                          merge_strategy="create_unique")
  counters = defaultdict(int)
  logs = list()
  for ftr in db.all_features():
    if len(list(db.parents(ftr.id))) == 0:
      if ftr.featuretype == "gene":
        _process_gene_feature(db, ftr, counters, logs)
      elif ftr.featuretype == "pseudogene":
        counters["pseudogene"] += 1
      elif ftr.featuretype == "region":
        counters["molecule"] += 1
      elif ftr.featuretype in ["sequence_feature"]:
        logs.append(("generic_sequence_feature", str(ftr)))
      elif ftr.featuretype == "riboswitch":
        counters["riboswitche"] += 1
      elif ftr.featuretype == "direct_repeat":
        if ftr.attributes["rpt_family"][0] == "CRISPR":
          counters["CRISPR"] += 1
        else:
          logs.append(f"direct_repeat_not_CRISPR\t{ftr}")
      else:
        logs.append(f"top_level_feature_other_type\t{ftr}")
  return _dict2list(counters, OUTPUT), logs

