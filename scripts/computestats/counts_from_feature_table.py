#!/usr/bin/env python3
"""
Count different types of features from a feature table
"""

count_functions = {
    "protein_coding_genes":
        lambda elems: elems[0] == "CDS" and elems[1] == "with_protein"
    }

