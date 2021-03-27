#!/usr/bin/env python3
import db_attributes
import db_create_tables
import db_load_results
from lib import db
from sqlalchemy import create_engine, select
from sqlalchemy.schema import MetaData
import yaml
from pathlib import Path
import uuid

ECHO=False
ENVVAR="PROSTTEST"
TESTDATA=Path("testdata")

def create_tables():
  args = db.args_from_env(ENVVAR)
  args["--verbose"] = ECHO
  for mod in ["computation_report", "plugin_description", "attribute"]:
    args["<file>"] = Path("dbschema")/f"{mod}.py"
    db_create_tables.main(args)

def run_create_attributes(definitions):
  args = db.args_from_env(ENVVAR)
  args["<definitions>"] = definitions
  args["--testmode"] = True
  args["--verbose"] = ECHO
  db_attributes.main(args)

def check_attributes(definitions):
  engine = create_engine(db.connstr_env(ENVVAR), echo=ECHO, future=True)
  with engine.connect() as connection:
    meta = MetaData()
    meta.reflect(bind=engine)
    adef_t = meta.tables["pr_attribute_definition"]
    rows = connection.execute(select(adef_t.c.name, adef_t.c.datatype,
                                     adef_t.c.computation_group)).all()
    assert(set(r.name for r in rows) == set(definitions.keys()))
    for r in rows:
      assert(r.datatype == definitions[r.name]["datatype"])
      assert(r.computation_group == definitions[r.name]["computation_group"])
    pfx = "pr_attribute_value_t"
    assert(f"{pfx}4" not in meta.tables)
    for n in [0, 1, 2, 3]:
      assert(f"{pfx}{n}" in meta.tables)
    assert(set(meta.tables[f"{pfx}0"].c.keys()) ==
           {"accession", "g1a_v", "g1a_c", "g1b_v", "g1b_c",
                         "g1c_v", "g1c_c", "g1_g"})
    assert(set(meta.tables[f"{pfx}1"].c.keys()) ==
           {"accession", "g1d_v", "g1d_c", "g1_g",
                         "g2a_v", "g2a_c", "g2b_v", "g2b_c", "g2_g"})
    assert(set(meta.tables[f"{pfx}2"].c.keys()) ==
           {"accession", "g3a_v0", "g3a_v1", "g3a_c", "g3_g"})
    assert(set(meta.tables[f"{pfx}3"].c.keys()) ==
           {"accession", "g3b_v0", "g3b_v1", "g3b_v2", "g3b_v3",
                         "g3b_c", "g3_g"})

def check_no_attributes(definitions):
  engine = create_engine(db.connstr_env(ENVVAR), echo=ECHO, future=True)
  with engine.connect() as connection:
    meta = MetaData()
    meta.reflect(bind=engine)
    adef_t = meta.tables["pr_attribute_definition"]
    rows = connection.execute(select(adef_t.c.name, adef_t.c.datatype,
                                     adef_t.c.computation_group)).all()
    assert(len(rows) == 0)
    for r in rows:
      assert(r.datatype == definitions[r.name]["datatype"])
      assert(r.computation_group == definitions[r.name]["computation_group"])
    pfx = "pr_attribute_value_t"
    assert(f"{pfx}4" not in meta.tables)
    for n in [0, 1, 2, 3]:
      assert(f"{pfx}{n}" in meta.tables)
      assert(set(meta.tables[f"{pfx}{n}"].c.keys()) == {"accession"})

def run_load_results(n):
  with open(TESTDATA/f"fake_report{n}.yaml") as f:
    args = db.args_from_env(ENVVAR)
    args["<report>"] = f
    args["<results>"] = TESTDATA/f"fake_results{n}.tsv"
    args["<plugin>"] = TESTDATA/f"fake_plugin{n}.py"
    db_load_results.main(args)

def run_destroy_attributes():
  args = db.args_from_env(ENVVAR)
  args["<definitions>"] = {}
  args["--drop"] = True
  args["--verbose"] = ECHO
  db_attributes.main(args)

def computation_id(n):
  return uuid.UUID(f'00000000-0000-0000-0000-00000000000{n}').bytes

def check_values_after_run_1():
  engine = create_engine(db.connstr_env(ENVVAR), echo=ECHO, future=True)
  with engine.connect() as connection:
    meta = MetaData()
    meta.reflect(bind=engine)
    t = [meta.tables[f"pr_attribute_value_t{n}"] for n in [0, 1, 2, 3]]
    assert(connection.execute(select(t[0])).all() ==
    [('A1', 1, None, computation_id(1), 2, None, 3, None),
     ('A2', 1, None, computation_id(1), 2, None, 3, None),
     ('A3', 1, None, computation_id(1), 2, None, 3, None),
     ('A4', 1, None, computation_id(1), 2, None, 3, None)])
    assert(connection.execute(select(t[1])).all() ==
    [('A1', 4, None, computation_id(1), None, None, None, None, None),
     ('A2', 4, None, computation_id(1), None, None, None, None, None),
     ('A3', 4, None, computation_id(1), None, None, None, None, None),
     ('A4', 4, None, computation_id(1), None, None, None, None, None)]
        )
    assert(connection.execute(select(t[2])).all() == [])
    assert(connection.execute(select(t[3])).all() == [])

create_tables()
with open(TESTDATA/"fake_attrs.yaml") as f:
  definitions = yaml.safe_load(f)
run_create_attributes(definitions)
check_attributes(definitions)
run_load_results(1)
check_values_after_run_1()
run_destroy_attributes()
check_no_attributes(definitions)
