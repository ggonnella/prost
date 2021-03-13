#!/usr/bin/env python3
def id_from_filename(filename):
  return "_".join(filename.split("/")[-1].split("_")[:2])

