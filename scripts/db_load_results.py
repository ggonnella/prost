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
             accession, attr_class, attr_instance, value[, score]
  report:    computation report file (yaml format)
  plugin:    plugin used for the computation

Options:
  --replace-plugin-record  replace existing db record for current version of the
                           plugin, if changed (default: fail if changed)
  --replace-report-record  replace existing db record for the computation
                           report, if changed (default: fail if changed)
  --verbose, -v    be verbose
  --version, -V    show script version
  --help, -h       show this help message
"""
from docopt import docopt
from schema import Schema, And, Or, Use, Optional
import os
from lib import snake, mod, plugins, db
import yaml
import MySQLdb
from sqlalchemy import create_engine, select, inspect
from sqlalchemy.orm import Session
from dbschema.attribute import AttributeDefinition, AttributeValueTables
from dbschema.plugin_description import PluginDescription
from dbschema.computation_report import ComputationReport

def different_fields(obj1, obj2):
  result = []
  for c in inspect(obj1).mapper.c:
    v1 = getattr(obj1, c.name)
    v2 = getattr(obj2, c.name)
    if v1 != v2:
      result.append((c, v1, v2))
  return result

def insert_update_or_compare(session, klass, newdata, primarykeys, replace, msg):
  newrow = klass(**newdata)
  oldrow = session.get(klass, primarykeys)
  if oldrow:
    if replace:
      session.merge(newrow)
    else:
      diff = different_fields(newrow, oldrow)
      if diff:
        raise RuntimeError(f"{msg}\nDifferences:\n"+
        "\n".join([f"  new {c}: {v1}\n  old {c}: {v2}" \
                   for c, v1, v2 in diff]))
  else:
    session.add(newrow)

def process_plugin_description(session, plugin, replace):
  insert_update_or_compare(session, PluginDescription, plugins.metadata(plugin),
      [plugin.ID, plugin.VERSION], replace,
      "Plugin metadata changed, without a version change\n"
      "Please use the --replace-plugin-record option or increase the "+\
      "plugin version number")

def check_plugin_key(data, plugin):
  for key in ["ID", "VERSION"]:
    report_value = data["plugin_"+key.lower()]
    exp_value = getattr(plugin, key)
    if report_value != exp_value:
      raise RuntimeError(f"Plugin {key} mismatch\n"+
          f"{key} from the plugin module: {exp_value}\n"+\
          f"{key} from the computation report: {report_value}\n")

def process_computation_report(session, reportfn, plugin, replace):
  data = yaml.safe_load(reportfn)
  check_plugin_key(data, plugin)
  insert_update_or_compare(session, ComputationReport, data, data["uuid"],
      replace, "The computation report was already stored in the database "+\
          "with a different content.\n"
          "Please either use the --replace-report-record option to "+\
          "replace the previous content.")
  return data["uuid"]

def main(args):
  engine = create_engine(db.connstr_from(args), echo=True, future=True)
  with engine.connect() as connection:
    with connection.begin():
      session = Session(bind=connection)
      plugin = mod.py_or_nim(args["<plugin>"], args["--verbose"])
      process_plugin_description(session, plugin,
                                 args["--replace-plugin-record"])
      computation_id = process_computation_report(session, args["<report>"],
                          plugin, args["--replace-report-record"])
      avt = AttributeValueTables(connection)
      avt.load_computation(computation_id, plugin.OUTPUT, args["<results>"])

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
                   Optional(str): object})
  return schema.validate(args)

if "snakemake" in globals():
  args = snake.args(snakemake,
        config=["<dbuser>", "<dbpass>", "<dbname>"],
        input=["<dbsocket>", "<results>", "<report>", "<plugin>"],
        params=["--replace-plugin-record", "--replace-report-record"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__, version="0.1")
  main(validated(args))
