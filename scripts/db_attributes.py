#!/usr/bin/env python3
"""
Create attribute definition records and attribute columns in the
attribute_value tables according to the definitions in a given YAML file.

Optionally remove attribute columns and definition records which are not
present in the YAML file.

Usage:
  db_attributes.py [options] <dbuser> <dbpass> <dbname> <dbsocket> <definitions>

Arguments:
  dbuser:       database user to use
  dbpass:       password of the database user
  dbname:       database name
  dbsocket:     connection socket file
  definitions:  YAML file containing the attribute definitions

Options:
  --drop           drop columns and delete definition record for all attributes
                   not present in the YAML file (be careful, DANGEROUS!)
  --check          check consistency of definition records and attribute columns
  --update         update definitions if changed
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
from docopt import docopt
from schema import Schema, And, Or, Use, Optional
from lib import snake
import os
import yaml
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker
from dbschema.attribute import AttributeDefinition, AttributeValueTables

def main(args):
  connstr = "".join(["mysql+mysqldb://", args["<dbuser>"],
                     ":", args["<dbpass>"], "@localhost/",
                     args["<dbname>"], "?unix_socket=",
                     args["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  avt = AttributeValueTables(engine)
  if args["--check"]:
    avt.check_consistency()
  if args["--update"]:
    Session = sessionmaker(engine, future=True)
    with Session.begin() as session:
      for adef in session.execute(select(AttributeDefinition)).scalars().all():
        if adef.name in args["<definitions>"]:
          for fieldname, fieldvalue in args["<definitions>"][adef.name].items():
            if fieldname == "datatype":
              if adef.datatype != fieldvalue:
                raise RuntimeError(f"Cannot update {adef.name} definition as "+\
                    "datatype changes are not allowed (previous value: "+\
                    f"{adef.datatype}, current value: {fieldvalue})")
            else:
              setattr(adef, fieldname, fieldvalue)
          session.add(adef)
  if args["--drop"]:
    to_delete = set(avt.attribute_names) - set(args["<definitions>"].keys())
    for aname in to_delete:
      avt.destroy_attribute(aname)
  for aname, adef in args["<definitions>"].items():
    if aname not in avt.attribute_names:
      adt = adef["datatype"]
      del adef["datatype"]
      avt.create_attribute(aname, adt, **adef)

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<definitions>": And(str, Use(open),
                                    Use(yaml.safe_load)),
                   "--update": Or(None, True, False),
                   "--drop": Or(None, True, False),
                   "--check": Or(None, True, False),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        config=["<dbuser>", "<dbpass>", "<dbname>"],
        input=["<dbsocket>", "<definitions>"],
        params=["--drop", "--check", "--update"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
