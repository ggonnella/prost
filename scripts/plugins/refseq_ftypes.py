#!/usr/bin/env python3
"""
Count features of different type and class
from a NCBI Refseq genome annotation GFF file.
"""

ID = "refseq_gff_count"
VERSION = "0.1.0"
INPUT = "genome annotation; from NCBI Refseq; gff, optionally gzipped"
OUTPUT = [
    "n_replicon",
    "n_protein_coding_gene",
    "n_tRNA",
    "n_rRNA_16S",
    "n_rRNA_5S",
    "n_rRNA_23S",
    "n_pseudogene",
    "n_riboswitch",
    "n_CRISPR",
    "n_translational_frameshift",
    "n_RNA_thermometer",
    "n_regulatory_promoter_element",
    "n_RNA_6S",
    "n_microsatellite",
    "n_monomeric_repeat",
    "n_tandem_repeat",
    "n_pseudoknot",
    "n_RNase_P_RNA",
    "n_SRP_RNA",
    "n_tmRNA",
    "n_ribosome_entry_site",
    "n_antisense_RNA",
    "n_hammerhead_ribozyme",
    "n_snoRNA"
]

METHOD = """
Count features based on terms defined in the Sequence ontology.

Replicons are counted as the number of features of type "region".
Genes are counted as "protein_coding_gene" if the contain at least
a child CDS feature.

Some feature types are assigned to more specific classes
(or prefixed by "other_" otherwise):
- rRNA  => rRNA_16S|5S|23S if product starts with "...S"
- ncRNA => RNA_6S          if product contains "6S RNA"
- mobile_genetic_element => insertion_sequence
    if mobile_element_type=insertion sequence:...
- direct_repeat => CRISPR if rpt_family=CRISPR
- binding_site => regulatory_promoter_element if Note contains: 'T_box_leader'
- repeat_region
    => CRISPR if Note=CRISPR repeat
    => microsatellite if rpt_family begins with SSR n, n > 1
    => monomeric_repeat if rpt_family begins with SSR n, n == 1
    => dispersed_repeat if Note contains: "Interspersed Repeat"
    => inverted_repeat if Note contains: "inverted"
    => other_direct_repeat if Note contains: "direct"
- sequence_feature
    => translational_frameshift if Note contains 'ribosomal frameshift element'
    => RNA_thermometer if Note contains: 'RNA thermometer'
    => regulatory_promoter_element if Note contains: 'leader region'
    => pseudoknot if Note contains: pseudoknot
Other features are counted by their feature_type column.

The counts for feature types contained in OUTPUT (n_<feature_type>)
are output in the results.

Other counts are output in the logs with the key "secondary_count" and are not
reliable or useful for comparisons. These include:
- features whose number is not interesting by itself: CDS, exon, intron
- very generic features: transcript
- generic features for which conditions for more specific assignment are
  described above, if those conditions are not met (prefixed in the output
  by "other_")
- features generally not annotated in Refseq (much less instances than expected
  are present in the annotations): origin_of_replication, pseudogenic_tRNA,
  operon, STS scRNA, recombination_feature, inverted_repeat, stem_loop,
  minus_10_signal, minus_35_signal, terminator, pseudogenic_tRNA,
  pseudogenic_rRNA, protein_binding_site
"""
IMPLEMENTATION = "Based on the py library Gffutils, using an in-memory database"
REQ_SOFTWARE = "pip library gffutils"

from collections import defaultdict
import gffutils


def _dict2list(d, keynames):
    [d.get(k, 0) for k in keynames]


def _process_gene_feature(db, gene, counters, logs):
    children = list(db.children(gene.id, level=1))
    if len(children) == 0:
        logs.append(f"unexpected\tgene_with_no_children\t{gene}")
        return
    child = children[0]
    if not all([child.featuretype == c.featuretype for c in children]):
        logs.append(f"unexpected\tgene_with_heterogeneous_children_type\t{gene}")
    if child.featuretype == "CDS":
        counters["n_protein_coding_gene"] += 1
    elif child.featuretype == "rRNA":
        product = child.attributes["product"][0]
        if product.startswith("16S"):
            counters["n_rRNA_16S"] += 1
        elif product.startswith("5S"):
            counters["n_rRNA_5S"] += 1
        elif product.startswith("23S"):
            counters["n_rRNA_23S"] += 1
        else:
            counters["n_other_rRNA"] += 1
    elif child.featuretype == "ncRNA":
        if "6S RNA" in child.attributes["product"][0]:
            counters["n_RNA_6S"] += 1
        else:
            counters["n_other_ncRNA"] += 1
    else:
        counters[f"n_{child.featuretype}"] += 1


def compute(filename, **kwargs):
    db = gffutils.create_db(filename, ":memory:", merge_strategy="create_unique")
    counters = defaultdict(int)
    logs = list()
    for ftr in db.all_features():
        if len(list(db.parents(ftr.id))) == 0:
            if ftr.featuretype == "region":
                counters["n_replicon"] += 1
            elif ftr.featuretype == "gene":
                _process_gene_feature(db, ftr, counters, logs)
            elif ftr.featuretype == "pseudogene":
                counters["n_pseudogene"] += 1
            elif ftr.featuretype in ["sequence_feature"]:
                if "ribosomal frameshift element" in ftr.attributes["Note"]:
                    counters["n_translational_frameshift"] += 1
                elif "RNA thermometer" in ftr.attributes["Note"]:
                    counters["n_RNA_thermometer"] += 1
                elif "leader region" in ftr.attributes["Note"]:
                    counters["n_regulatory_promoter_element"] += 1
                elif "pseudoknot" in ftr.attributes["Note"]:
                    counters["n_pseudoknot"] += 1
                else:
                    counters["n_other_sequence_feature"] += 1
            elif ftr.featuretype == "riboswitch":
                counters["n_riboswitch"] += 1
            elif ftr.featuretype == "mobile_genetic_element":
                if ftr.attributes["mobile_element_type"][0].startswith(
                    "insertion sequence"
                ):
                    counters["n_insertion_sequence"] += 1
                else:
                    counters["n_other_mobile_genetic_element"] += 1
            elif ftr.featuretype == "direct_repeat":
                if ftr.attributes["rpt_family"][0] == "CRISPR":
                    counters["n_CRISPR"] += 1
                else:
                    counters["n_other_direct_repeat"] += 1
            elif ftr.featuretype == "repeat_region":
                if ftr.attributes["rpt_family"][0] == "SSR 1mer":
                    counters["n_monomeric_repeat"] += 1
                elif ftr.attributes["rpt_family"][0].startswith("SSR"):
                    counters["n_microsatellite"] += 1
                elif "CRISPR" in ftr.attributes["Note"][0]:
                    counters["n_CRISPR"] += 1
                elif "direct" in ftr.attributes["Note"][0]:
                    counters["n_other_direct_repeat"] += 1
                elif "inverted" in ftr.attributes["Note"][0]:
                    counters["n_inverted_repeat"] += 1
                elif "Interespersed Repeat" in ftr.attributes["Note"][0]:
                    counters["n_dispersed_repeat"] += 1
                else:
                    counters["n_other_repeat_region"] += 1
            elif ftr.featuretype == "binding_site":
                if "T-box leader" in ftr.attributes["Note"]:
                    counters["n_regulatory_promoter_element"] += 1
                else:
                    counters["n_other_binding_site"] += 1
            else:
                counters[f"n_{ftr.feature_type}"] += 1
        for k, v in counters.items():
            if k not in OUTPUT:
                logs.append(f"secondary_count\t{k}\t{v}")
    return _dict2list(counters, OUTPUT), logs
