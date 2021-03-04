#!/usr/bin/env python3
"""
Organize the content of a collection of JSON objects
in a tabular way, so that it can be loaded into a relational database.

Requirements:
  The input file is a TSV, in which one column contains the Json
  and, optionally, another column contains further data to include
  in all tables (e.g. an ID).

  The JSON must have a structure like the following:
    dict:
     dict:
      ... (any number of dict nesting levels)
       list:
         dict:
           scalar values
  A table is generated from each "list" level.
  All dict keys must match [a-zA-Z]+[0-9a-zA-Z_]*.

Usage:
  json_to_tabular.py [options] <tsv> <colnum> <pfx> [<idcolnum> <idcol>]

Arguments:
  tsv:      tsv file, with a column containing JSON data
  colnum:   1-based number of column containing JSON data
  pfx:      prefix for the table names

Optional arguments:
  idcolnum: 1-based number of column containing IDs
  idcol:    column name to use for the IDs (do not use "id")

Output is a TSV file with lines containing:
  (1) table name
  (2) column name
  (3) value (JSON format)
Rows are divided by lines containing only '###'

Options:
  --verbose, -v  be verbose
  --version, -V  show script version
  --help, -h     show this help message
"""

from docopt import docopt
from schema import Schema, Use, And, Or
import json

NormalizedColnames = {
    "class": "klass"
  }

def handle_list(tablename, rows, id_data, idcol):
  for row in rows:
    if id_data is not None:
      print("\t".join([tablename, idcol, id_data]))
    for colname, value in row.items():
      print("\t".join([tablename,
        NormalizedColnames.get(colname, colname),
        json.dumps(value)]))
    print("###")

def handle_dict(pfx, json_data, id_data, idcol):
  for k, v in json_data.items():
    if isinstance(v, dict):
      shortk = "".join([w[0] for w in k.split("_")])
      pfx2 = pfx + "_" + shortk
      handle_dict(pfx2, v, id_data, idcol)
    else:
      assert(isinstance(v, list))
      tablename = pfx + "_" + k
      handle_list(tablename, v, id_data, idcol)

def main(arguments):
  for line in arguments["<tsv>"]:
    elems = line.rstrip().split("\t")
    if arguments["<idcolnum>"] is not None:
      id_data = elems[arguments["<idcolnum>"]]
    else:
      id_data = None
    json_data = json.loads(elems[arguments["<colnum>"]])
    handle_dict(arguments["<pfx>"], json_data, id_data, arguments["<idcol>"])

def validated(arguments):
  colnum = And(Use(int), lambda i: i>0, Use(lambda i: i-1))
  schema = Schema({"<tsv>": Use(open),
                  "<colnum>": colnum,
                  "<pfx>": len,
                  "<idcol>": Or(None, len),
                  "<idcolnum>": Or(None, colnum)},
                  ignore_extra_keys=True)
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(validated(arguments))
