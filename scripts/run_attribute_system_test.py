#!/usr/bin/env python3
import db_attributes
import db_create_tables
import db_load_results
from lib import db
from sqlalchemy import create_engine, select, text
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
    args["--verbose"] = ECHO
    db_load_results.main(args)

def run_destroy_attributes():
  args = db.args_from_env(ENVVAR)
  args["<definitions>"] = {}
  args["--drop"] = True
  args["--verbose"] = ECHO
  db_attributes.main(args)

def computation_id(n):
  return uuid.UUID(f'00000000-0000-0000-0000-00000000000{n}').bytes

def get_attribute_value_rows():
 engine = create_engine(db.connstr_env(ENVVAR), echo=ECHO, future=True)
 with engine.connect() as connection:
   meta = MetaData()
   meta.reflect(bind=engine)
   tnums = [0,1,2,3]
   r = []
   for n in tnums:
     t = meta.tables[f"pr_attribute_value_t{n}"]
     rows = connection.execute(select(t)).all()
     r.append({row.accession: row for row in rows})
   return r

def check_g1_t0(row, run):
  if run != 2:
    assert(row.g1a_v == 1)
    assert(row.g1b_v == 2)
    assert(row.g1c_v == 3)
    # computation ID in _g and not _c, as whole grp computed
    assert(row.g1_g == computation_id(1))
    for field in ["g1a_c", "g1b_c", "g1c_c"]:
      assert(getattr(row, field) == None)
  else:
    # computation IDs stored in _c for re-computed values
    # as not entire group recomputed
    assert(row.g1a_v == 5)
    assert(row.g1a_c == computation_id(2))
    assert(row.g1b_v == 2)
    assert(row.g1b_c == None)
    assert(row.g1_g == computation_id(1))
    assert(row.g1c_v == 6)
    assert(row.g1c_c == computation_id(2))

def check_g1_t1(row, run):
  assert(row.g1d_c == None)
  if run != 2:
    assert(row.g1d_v == 4)
    assert(row.g1_g == computation_id(1))
  else:
    assert(row.g1d_v == 7)
    assert(row.g1_g == computation_id(2))

def check_g2_t1(row, run):
  if run < 3:
    for field in ["g2a_v", "g2a_c", "g2b_v", "g2b_c", "g2_g"]:
      assert(getattr(row, field) == None)
  else:
    assert(row.g2a_v == 8)
    assert(row.g2a_c == None)
    assert(row.g2_g == computation_id(3))
    if run == 3:
      assert(row.g2b_v == 9)
      assert(row.g2b_c == None)
    else:
      assert(row.g2b_v == 10)
      assert(row.g2b_c == computation_id(4))

def check_g3_t2(row, run):
  assert(row.g3a_v0 == 10)
  assert(row.g3a_v1 == 2.3)
  assert(row.g3a_c == None)
  if run <= 6:
    assert(row.g3_g == computation_id(5))
  else:
    assert(row.g3_g == computation_id(7))

def check_g3_t3(row, run):
  assert(row.g3b_v0 == 1.0)
  assert(row.g3b_v1 == 11)
  assert(row.g3b_v2 == 12)
  assert(row.g3b_v3 == 13)
  assert(row.g3b_c == None)
  if run == 5:
    assert(row.g3_g == computation_id(5))
  else:
    assert(row.g3_g == computation_id(6))

def check_values_after_run(run):
    r = get_attribute_value_rows()
    for accession in ['A1','A2','A3','A4']:
      check_g1_t0(r[0][accession], run)
      check_g1_t1(r[1][accession], run)
      check_g2_t1(r[1][accession], run)
    if run <= 4:
      assert(r[2] == {})
      assert(r[3] == {})
    else:
      for accession in ['A1','A2','A3','A4']:
        check_g3_t2(r[2][accession], run)
        check_g3_t3(r[3][accession], run)

def drop_dangling_attribute_tables():
  engine = create_engine(db.connstr_env(ENVVAR), echo=ECHO, future=True)
  with engine.connect() as connection:
    connection.execute(text("DROP TABLE IF EXISTS pr_attribute_definition"))
    for n in [0,1,2,3, "temporary"]:
      connection.execute(text(f"DROP TABLE IF EXISTS pr_attribute_value_t{n}"))

drop_dangling_attribute_tables()
create_tables()
with open(TESTDATA/"fake_attrs.yaml") as f:
  definitions = yaml.safe_load(f)
run_create_attributes(definitions)
check_attributes(definitions)
run_load_results(1)
check_values_after_run(1)
run_load_results(2)
check_values_after_run(2)
run_load_results(1)
check_values_after_run(1)
run_load_results(3)
check_values_after_run(3)
run_load_results(4)
check_values_after_run(4)
run_load_results(5)
check_values_after_run(5)
run_load_results(6)
check_values_after_run(6)
run_load_results(7)
check_values_after_run(7)
run_destroy_attributes()
check_no_attributes(definitions)
