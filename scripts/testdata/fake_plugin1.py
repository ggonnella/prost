#!/usr/bin/env python3
"""
Does nothing, for test purposes
"""

ID =      "fake1"
VERSION = "0.1.0"
INPUT = "test"
OUTPUT =  ["g1a", "g1b", "g1c", "g1d"]

def compute(filename, **kwargs):
  return [1,2,3,4], ["test\tthis is a test by fake1"]

def fake_data():
  for acc in ["A1","A2","A3","A4"]:
    print("\t".join([acc]+[str(e) for e in compute(acc)[0]]))

if __name__ == "__main__":
  fake_data()