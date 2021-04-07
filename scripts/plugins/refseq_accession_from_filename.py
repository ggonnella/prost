#!/usr/bin/env python3
def compute_id(filename):
  return "_".join(filename.split("/")[-1].split("_")[:2])

