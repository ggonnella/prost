#!/usr/bin/env python3
"""
Vizualize information about the structure of a set of JSON records
in a TSV file.

Usage:
  ./summarize_json_records.py [options] <tsvfile> <colnum>

Arguments:
  tsvfile  filename containing the JSON information
  colnum   1-based column number to use

Options:
  --verbose, -v   be verbose
  --help, -h      show this help message
  --version, -V   show the script version
"""

from docopt import docopt
from schema import Schema, Use, And
from collections import Counter
import json
import sys

def json_collect_info(stdata, stinfo, key):
  if key not in stinfo:
    stinfo[key] = dict()
  if "#" not in stinfo[key]:
    stinfo[key]["#"] = 0
  stinfo[key]["#"] += 1
  if isinstance(stdata, dict):
    for k, v in stdata.items():
      json_collect_info(v, stinfo[key], k)
  elif isinstance(stdata, list):
    if "[]" not in stinfo[key]:
      stinfo[key]["[]"] = {"__len__": []}
    stinfo[key]["[]"]["__len__"].append(len(stdata))
    for e in stdata:
      json_collect_info(e, stinfo[key], "[]")
  else:
    if stdata not in stinfo[key]:
      stinfo[key][stdata] = 0
    stinfo[key][stdata] += 1

def format_key(k, maxklen):
  k = str(k)
  klen=len(k)
  k = k.replace("\r", " ")
  k = k.replace("\n", " ")
  if len(k) > maxklen:
    k = k[:(maxklen-3)]+f"... (len:{klen})"
  return k

def print_info(stinfo, level):
  pfx = "  "*level
  maxklen=40
  for k, v in stinfo.items():
    if k == "#":
      continue
    elif k == "[]":
      print(f"{pfx}<list>:")
      sumlen = sum(v["__len__"])
      avglen = sumlen / len(v["__len__"])
      minlen = min(v["__len__"])
      nminlen = v["__len__"].count(minlen)
      maxlen = max(v["__len__"])
      if minlen == 0 and maxlen > 0:
        l2 = [x for x in v["__len__"] if x > 0]
        minlen2 = min(l2)
        nminlen2 = l2.count(minlen2)
      nmaxlen = v["__len__"].count(maxlen)
      if minlen == maxlen:
        print(f"{pfx}  <len>: {minlen} (constant)")
      else:
        if minlen == 0:
          if minlen2 == maxlen:
            print(f"{pfx}  <len>: 0(n={nminlen}) or {minlen2}(n=nminlen2)")
          else:
            print(f"{pfx}  <len>: 0(n={nminlen}) or "+\
                f"{minlen2}(n={nminlen2}) .. {maxlen}(n={nmaxlen}) "+\
                f"(avg: {avglen:.2f}, total:{sumlen})")
        else:
          print(f"{pfx}  <len>: {minlen}(n={nminlen}) .. {maxlen}(n={nmaxlen}) "+\
              f"(avg: {avglen:.2f}, total:{sumlen})")
      print(f"{pfx}  <elements>:")
      del v["__len__"]
      print_info(v, level+2)
    elif isinstance(v, dict):
      print(f"{pfx}{format_key(k, maxklen)}: {v['#']}")
      print_info(v, level+1)
  leaves={k:v for (k, v) in stinfo.items() if isinstance(v, int)}
  if "#" in leaves:
    del leaves["#"]
  leaves = Counter(leaves)
  maxstinfolen=10
  if len(leaves) <= maxstinfolen:
    for k, v in leaves.most_common():
      print(f"{pfx}{format_key(k, maxklen)}: {v}")
  else:
    remaining = sum(leaves.values())
    for k, v in leaves.most_common(maxstinfolen//2):
      print(f"{pfx}{format_key(k, maxklen)}: {v}")
      remaining -= v
    rdiff = len(leaves) - (maxstinfolen//2)
    print(f"{pfx}... and other {rdiff} different values "+\
        f"in {remaining} elements")

def main(arguments):
  info = dict()
  for line in arguments["<tsvfile>"]:
    elems = line.rstrip().split("\t")
    jdata = json.loads(elems[arguments["<colnum>"]])
    json_collect_info(jdata, info, "\\")
  print_info(info, 0)

def validated(arguments):
  colnum = And(Use(lambda i: int(i)-1), lambda i: i>=0)
  schema = Schema({"<tsvfile>": Use(open), "<colnum>": colnum},
                  ignore_extra_keys=True)
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(validated(arguments))
