#!/usr/bin/env python3
"""
Create attribute definition records and attribute columns in the
attribute_value tables according to the definitions in a given YAML file.

Optionally remove attribute columns and definition records which are not
present in the YAML file.

Usage:
  db_attributes.py [options] {db_args_usage} <definitions>

Arguments:
{db_args}
  definitions:  YAML file containing the attribute definitions

Options:
  --drop           drop columns and delete definition record for all attributes
                   not present in the YAML file (be careful, DANGEROUS!)
  --check          check consistency of definition records and attribute columns
  --update         update definitions if changed
  --testmode       use the parameters for tests
{common}
"""
from schema import And, Or, Use
from lib import db, scripts
import yaml
from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session
from dbschema.attribute import AttributeDefinition, AttributeValueTables
import snacli

def update(connection, definitions):
  if definitions:
    session = Session(connection)
    select_adefs = select(AttributeDefinition)
    for adef in session.execute(select_adefs).scalars().all():
      if adef.name in definitions:
        for fname, fvalue in definitions[adef.name].items():
          if fname == "datatype":
            if adef.datatype != fvalue:
              raise RuntimeError(f"Cannot update {adef.name} definition "+\
                  "as datatype changes are not allowed (previous: "+\
                  f"{adef.datatype}, current: {fvalue})")
          else:
            setattr(adef, fname, fvalue)
        session.add(adef)
    session.commit()

def drop(avt, definitions):
  for aname in avt.attribute_names:
    if not definitions or aname not in definitions:
      avt.destroy_attribute(aname)

def insert(avt, definitions):
  if definitions:
    for aname, adef in definitions.items():
      if aname not in avt.attribute_names:
        adef = adef.copy()
        adt = adef["datatype"]
        del adef["datatype"]
        avt.create_attribute(aname, adt, **adef)

def main(args):
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"],
                         future=True)
  with engine.connect() as connection:
    with connection.begin():
      avt = AttributeValueTables(connection)
      if args["--testmode"]: avt.TARGET_N_COLUMNS = 9
      if args["--check"]:  avt.check_consistency()
      if args["--update"]: update(connection, args["<definitions>"])
      if args["--drop"]:   drop(avt, args["<definitions>"])
      insert(avt, args["<definitions>"])

def validated(args):
  return scripts.validate(args, db.args_schema,
      {"<definitions>": And(str, Use(open), Use(yaml.safe_load)),
       "--update": Or(None, True, False), "--testmode": Or(None, True, False),
       "--drop": Or(None, True, False), "--check": Or(None, True, False)})

with snacli.args(db.snake_args, input=["<definitions>"],
                    params=["--drop", "--check", "--update", "--testmode",
                            "--verbose"],
                 docvars={"common": scripts.args_doc,
                 "db_args": db.args_doc, "db_args_usage": db.args_usage},
                 version="0.1") as args:
  if args: main(validated(args))
