#!/usr/bin/env python3
"""
Count different types of features from a feature table
"""

counters = {
    "protein_coding_genes":
        lambda r: r["feature"] == "CDS" and \
                  r["class"]   == "with_protein"
    }

