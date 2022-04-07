#!/usr/bin/env python3
"""
Load computation results into DB from a TSV file and a computation report.

Usage:
  db_load_results.py [options] {db_args_usage}
                               <results> <report> <plugin>

Arguments:
{db_args}
  results:   results file, tsv with columns:
             accession, attr_class, attr_instance, value[, score]
  report:    computation report file (yaml format)
  plugin:    plugin used for the computation

Options:
  --replace-plugin-record  replace existing db record for current version of the
                           plugin, if changed (default: fail if changed)
  --replace-report-record  replace existing db record for the computation
                           report, if changed (default: fail if changed)
{common}
"""
from docopt import docopt
from schema import And, Or, Use
from lib import snake, plugins, db, scripts
import multiplug
import yaml
import os
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import Session
from dbschema.attribute import AttributeValueTables
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

def insert_update_or_compare(session, klass, newdata,
                             primarykeys, replace, msg):
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
  insert_update_or_compare(session, PluginDescription,
      plugins.plugin_metadata_str(plugin),
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
  if os.stat(args["<results>"]).st_size == 0:
    if args["--verbose"]:
      sys.stderr.write("# nothing to load, as results file is empty \n")
    return
  engine = create_engine(db.connstr_from(args), echo=args["--verbose"],
                         future=True)
  with engine.connect() as connection:
    with connection.begin():
      session = Session(bind=connection)
      plugin = multiplug.importer(args["<plugin>"], verbose=args["--verbose"],
                                  **plugins.COMPUTE_PLUGIN_INTERFACE)
      process_plugin_description(session, plugin,
                                 args["--replace-plugin-record"])
      computation_id = process_computation_report(session, args["<report>"],
                          plugin, args["--replace-report-record"])
      session.commit()
      avt = AttributeValueTables(connection)
      avt.load_computation(computation_id, plugin.OUTPUT, args["<results>"])

def validated(args):
  return scripts.validate(args, db.args_schema,
                  {"<results>": And(str, open),
                   "<report>": And(str, Use(open)),
                   "<plugin>": And(str, open),
                   "--replace-plugin-record": Or(None, True, False),
                   "--replace-report-record": Or(None, True, False)})

if "snakemake" in globals():
  args = snake.args(snakemake, db.snake_args,
        input=["<results>", "<report>", "<plugin>"],
        params=["--replace-plugin-record", "--replace-report-record",
                "--verbose"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(db_args=db.args_doc, db_args_usage=db.args_usage,
               common=scripts.args_doc), version="0.1")
  main(validated(args))
