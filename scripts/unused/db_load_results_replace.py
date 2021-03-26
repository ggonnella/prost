#!/usr/bin/env python3
"""
Load computation results into DB from a TSV file and a computation report.

Usage:
  db_load_results.py [options] <dbuser> <dbpass> <dbname> <dbsocket>
                               <results> <report> <plugin>

THIS APPROACH DOES NOT WORK, AS IT REPLACES THE ENTIRE ROW DESTROYING
EXISTING DATA!

FURTHERMORE IT MIXES MYSQLDB AND SQLALCHEMY

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  results:   results file, tsv with columns:
             accession, attr_class, attr_instance, value[, score]
  report:    computation report file (yaml format)
  plugin:    plugin used for the computation

Options:
  --replace-plugin-record  replace existing db record for current version of the
                           plugin, if changed (default: fail if changed)
  --replace-report-record  replace existing db record for the computation
                           report, if changed (default: fail if changed)
  --scored, -s     results file has a float score column (default: False)
  --ignore, -i     use IGNORE on repeated primary key (default: REPLACE)
  --dropkeys, -k   drop non-unique indices before inserting
                   and re-compute them after inserting
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
from docopt import docopt
from schema import Schema, And, Or, Use, Optional
import os
from lib import snake, mysql, sqlwriter, mod, plugins
import yaml
import MySQLdb
from sqlalchemy import create_engine, select
from sqlalchemy.orm import sessionmaker
from dbschema.attribute import AttributeDefinition, AttributeValueTables

def insert_sql(tablename, columns, replace_columns=None):
  statement = f"INSERT INTO {tablename} ("
  statement += ", ".join(columns)
  statement += ") VALUES("
  statement += ", ".join(["%("+c+")s" for c in columns])
  statement += ")"
  if replace_columns:
    statement +=" ON DUPLICATE KEY UPDATE "
    statement += ", ".join(f"{c}=VALUES({c})" for c in replace_columns)
  statement += ";"
  return statement

def check_data_consistency(existing, data, desc, advice):
  for k, v in existing.items():
    if v != data.get(k, v):
      raise RuntimeError(f"Inconsistency in {desc} data\n"+\
          f"Field '{k}' in existing record: {v}\n"+\
          f"Field '{k}' in provided file: {data[k]}\n"+\
          advice)

def insert_or_check(cursor, tablename, data, primary_key, replace,
                    desc, advice):
  columns = list(data.keys())
  if replace:
    replace = columns.copy()
    for key in primary_key:
      replace.remove(key)
  statement = insert_sql(tablename, columns, replace)
  try:
    cursor.execute(statement, data)
  except MySQLdb.IntegrityError:
    where = " AND ".join([f"{k} = %({k})s" for k in primary_key])
    cursor.execute(f"SELECT * from {tablename} WHERE {where};", data)
    check_data_consistency(cursor.fetchone(), data, desc, advice)
  return data

def process_plugin_description(cursor, plugin, replace):
  insert_or_check(cursor, "pr_plugin_description",
                  plugins.metadata(plugin), ["ID", "VERSION"],
                  replace, "plugin description",
                  "Please either use the --replace-plugin-record option or "+\
                  "increase the plugin version number\n")

def check_plugin_key_consistency(data, plugin):
  for key in ["ID", "VERSION"]:
    report_value = data["plugin_"+key.lower()]
    exp_value = getattr(plugin, key)
    if report_value != exp_value:
      raise RuntimeError(f"Plugin {key} mismatch\n"+
          f"{key} from the plugin module: {exp_value}\n"+\
          f"{key} from the computation report: {report_value}\n")

def process_computation_report(cursor, reportfn, exp_plugin_data, replace):
  data = yaml.safe_load(reportfn)
  check_plugin_key_consistency(data, exp_plugin_data)
  insert_or_check(cursor, "pr_computation_report", data, ["uuid"], replace,
                  "computation report",
                  "Please either use the --replace-report-record option\n")
  return data["uuid"]

def locate_attr_columns(args, attribute_names):
  connstr = "".join(["mysql+mysqldb://", args["<dbuser>"],
                     ":", args["<dbpass>"], "@localhost/",
                     args["<dbname>"], "?unix_socket=",
                     args["<dbsocket>"]])
  engine = create_engine(connstr, echo=True)
  avt = AttributeValueTables(engine)
  return avt.locations_for_attributes(attribute_names)

def main(args):
  db = mysql.connect(args)
  cursor = db.cursor(MySQLdb.cursors.DictCursor)
  plugin = mod.importer(args["<plugin>"], args["--verbose"])
  process_plugin_description(cursor, plugin,
                   args["--replace-plugin-record"])
  computation_id = process_computation_report(cursor, args["<report>"],
                      plugin, args["--replace-report-record"])
  locations = locate_attr_columns(args, plugin.OUTPUT)
  for tablename, tlocations in locations["tables"].items():
    fixed_data = {cn: computation_id for cn in tlocations["ccols_to_set"]}
    for cn in tlocations["ccols_to_unset"]: fixed_data[cn] = None
    skipfields = [i+1 for i, cn in enumerate(locations["vcols"]) \
                  if not cn in tlocations["vcols_to_set"]]
    statements = sqlwriter.load_data_sql(args["<results>"], tablename,
                      ["accession"]+tlocations["vcols_to_set"],
                      skipfields, fixed_data, args["--ignore"],
                      args["--dropkeys"], False, "")
    for statement in statements:
      print(statement)
      #cursor.execute(statement, fixed_data)
  db.commit()
  cursor.close()

def validated(args):
  schema = Schema({"<dbuser>": And(str, len),
                   "<dbpass>": And(str, len),
                   "<dbname>": And(str, len),
                   "<dbsocket>": And(str, len, os.path.exists),
                   "<results>": And(str, open),
                   "<report>": And(str, Use(open)),
                   "<plugin>": And(str, open),
                   "--replace-plugin-record": Or(None, True, False),
                   "--replace-report-record": Or(None, True, False),
                   "--ignore": Or(None, True, False),
                   "--dropkeys": Or(None, True, False),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        config=["<dbuser>", "<dbpass>", "<dbname>"],
        input=["<dbsocket>", "<results>", "<report>", "<plugin>"],
        params=["--replace-plugin-record", "--replace-report-record",
                "--ignore", "--dropkeys"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
