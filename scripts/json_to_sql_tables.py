#!/usr/bin/env python3
"""
Generate code for SqlAlchemy tables to contain the data
from a list of JSON objects.

Limitation:
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
  ./json_to_tables.py [options] <tsv> <colnum> <pfx> <idcol>

Arguments:
  tsv:      tsv file, with a column containing JSON data
  colnum:   1-based number of column containing JSON data
  pfx:      prefix for the table names
  idcol:    name of the column containing IDs

Options:
  --verbose, -v  be verbose
  --version, -V  show script version
  --help, -h     show this help message
"""

from docopt import docopt
from schema import Schema, Use, And
import json
import textwrap

f = open("/home/gonnella/data/bacdive/bacdive.details")
jsoncol = 3

def emit_prelude(arguments):
  print(textwrap.dedent("""\
  #!/usr/bin/env python3
  \"\"\"
  This file was autogenerated, based on the content of file:\
  """))
  print(arguments["<tsv>"])
  print(textwrap.dedent("""\
  \"\"\"

  from sqlalchemy.sql import func
  from sqlalchemy.ext.declarative import declarative_base
  from sqlalchemy import Column, Integer, String, Sequence, \\
                         DateTime, Text, Boolean, Float
  from sqlalchemy_repr import PrettyRepresentableBase

  Base = declarative_base(cls=PrettyRepresentableBase)

  utf8_cs_args = {'mysql_charset': 'utf8', 'mysql_collate': 'utf8_bin'}
  """))

def normalize_colname(colname):
  if colname == "class":
    return "klass"
  return colname

def handle_list(data, tablename, rows):
  for row in rows:
    for colname, value in row.items():
      colname = normalize_colname(colname)
      if colname not in data[tablename]:
        data[tablename][colname] = []
      if value is not None:
        data[tablename][colname].append(value)

def handle_dict(data, pfx, d):
  for k, v in d.items():
    if isinstance(v, dict):
      shortk = "".join([w[0] for w in k.split("_")])
      pfx += "_" + shortk
      handle_dict(data, pfx, v)
    else:
      assert(isinstance(v, list))
      tablename = pfx + "_" + k
      if tablename not in data:
        data[tablename] = {}
      handle_list(data, tablename, v)

def parse_file(arguments):
  result = {}
  with open(arguments["<tsv>"]) as f:
    for line in f:
      elems = line.rstrip().split("\t")
      json_data = json.loads(elems[arguments["<colnum>"]])
      handle_dict(result, arguments["<pfx>"], json_data)
  return result

def emit_tables_code(data, arguments):
  for tablename in data.keys():
    klassname = "".join(w.title() for w in tablename.split("_"))
    print("")
    print(f"class {klassname}(Base):")
    print(f"  __tablename__ = '{tablename}'")
    print("  id = Column(Integer, Sequence('{tablename}_id_seq'),"+
          " primary_key=True)")
    print(f"  {arguments['<idcol>']} = Column(Integer)")
    print("  time_updated = Column(DateTime,")
    print("    server_default=func.now(), onupdate=func.now())")
    for colname, coldata in data[tablename].items():
      if all(isinstance(v, bool) for v in coldata):
        print(f"  {colname} = Column(Boolean)")
      elif all(isinstance(v, int) for v in coldata):
        print(f"  {colname} = Column(Integer)")
      elif all(isinstance(v, float) for v in coldata):
        print(f"  {colname} = Column(Float)")
      else:
        lens = [len(str(v)) for v in coldata]
        maxlen = max(lens)
        if maxlen < 32:
          print(f"  {colname} = Column(String(64))")
        elif maxlen < 64:
          print(f"  {colname} = Column(String(128))")
        elif maxlen < 128:
          print(f"  {colname} = Column(String(256))")
        else:
          print(f"  {colname} = Column(Text({maxlen*2}))")
    print("  __table_args__ = utf8_cs_args")

def main(arguments):
  emit_prelude(arguments)
  data = parse_file(arguments)
  emit_tables_code(data, arguments)

def validated(arguments):
  schema = Schema({"<tsv>": open,
                  "<colnum>": And(Use(int), lambda i: i>0, Use(lambda i: i-1)),
                  "<pfx>": len,
                  "<idcol>": len}, ignore_extra_keys=True)
  return schema.validate(arguments)

if __name__ == "__main__":
  arguments = docopt(__doc__, version="0.1")
  main(validated(arguments))
