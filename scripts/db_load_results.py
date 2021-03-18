#!/usr/bin/env python3
"""
Load computation results into DB from a TSV file and a computation report.

Usage:
  db_load_results.py [options] <dbuser> <dbpass> <dbname> <dbsocket>
                               <results> <report> <plugin>

Arguments:
  dbuser:    database user to use
  dbpass:    password of the database user
  dbname:    database name
  dbsocket:  connection socket file
  results:   results file, tsv with columns:
             accession, attr_class, attr_instance, value
  report:    computation report file (yaml format)
  plugin:    plugin used for the computation

Options:
  --replace-plugin-record  replace the db record for current version of the
                           plugin, if changed (default: fail if changed)
  --float, -f      value is float (default: int)
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
from lib import snake, mysql, sqlwriter, mod
import yaml
import MySQLdb

def process_plugin_description(cursor, plugin, replace):
  yamldata = plugin.__doc__.split("\n---\n")[1]
  data = yaml.safe_load(yamldata)
  columns = list(data.keys())
  statement = "INSERT INTO plugin("
  statement += ", ".join(columns)
  statement += ") VALUES("
  statement += ", ".join(["%("+c+")s" for c in columns])
  statement += ")"
  if replace:
    non_primary = columns.copy()
    non_primary.remove("id")
    non_primary.remove("version")
    statement +=" ON DUPLICATE KEY UPDATE "
    statement += ", ".join(f"{c}=VALUES({c})" for c in non_primary)
  try:
    cursor.execute(statement, data)
  except MySQLdb.IntegrityError:
    cursor.execute("SELECT * from plugin WHERE id = %(id)s AND "+\
        "version = %(version)s", data)
    existing = cursor.fetchone()
    for k, v in existing.items():
      if v != data.get(k, None):
        raise RuntimeError("Plugin description changed\n"+\
            f"Field '{k}' in plugin description: {data[k]}\n"+\
            f"Field '{k}' in existing record: {data[k]}\n"+\
            "Please either use the --replace-plugin-record option or "+\
            "increase the plugin version number\n")
  return data

def process_computation_report(cursor, reportfn, exp_plugin_data):
  data = yaml.safe_load(reportfn)
  for key in ["id", "version"]:
    report_value = data["plugin_"+key]
    exp_value = exp_plugin_data[key]
    if report_value != exp_value:
      raise RuntimeError(f"Plugin {key} mismatch\n"+
          f"{key} from the plugin module: {exp_value}\n"+\
          f"{key} from the computation report: {report_value}\n")
  # insert data
  columns = data.keys()
  statement = "INSERT INTO computation("
  statement += ", ".join(columns)
  statement += ") VALUES("
  statement += ", ".join(["%("+c+")s" for c in columns])
  statement += ");"
  cursor.execute(statement, data)
  cursor.execute("SELECT LAST_INSERT_ID() AS ID;")
  computation_id = cursor.fetchone()["ID"]
  return computation_id

def main(args):
  tablename = "assembly_float_value" if args["--float"] else \
              "assembly_int_value"
  columns = ["accession", "attr_class", "attr_instance", "value"]
  db = mysql.connect(args)
  cursor = db.cursor(MySQLdb.cursors.DictCursor)
  plugin = mod.importer(args["<plugin>"], args["--verbose"])
  plugin_data = \
    process_plugin_description(cursor, plugin, args["--replace-plugin-record"])
  computation_id = \
    process_computation_report(cursor, args["<report>"], plugin_data)
  fixed_data = {"computation_id": computation_id}
  statements = sqlwriter.load_data_sql(args["<results>"], tablename, columns,
                    fixed_data, args["--ignore"], args["--dropkeys"], False, "")
  for statement in statements:
    cursor.execute(statement)
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
                   "--float": Or(None, True, False),
                   "--ignore": Or(None, True, False),
                   "--dropkeys": Or(None, True, False),
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        config=["<dbuser>", "<dbpass>", "<dbname>"],
        input=["<dbsocket>", "<results>", "<report>", "<plugin>"],
        params=["--replace-plugin-record", "--float", "--ignore", "--dropkeys"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
