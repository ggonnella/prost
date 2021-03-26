"""
Helper methods for working with the database
"""

from sqlalchemy.engine.url import URL
import os

DB_DRIVER="mysql+mysqldb"
DB_HOST="localhost"

def _connstr(u,p,d,s):
  if not(s):
    raise RuntimeError(f"DB unix socket not provided")
  elif not os.path.exists(s):
    raise RuntimeError(f"DB unix socket does not exist: {s}")
  elif not(u):
    raise RuntimeError("Database user name not provided")
  elif not(p):
    raise RuntimeError(f"Database password for user '{u}' not provided")
  elif not(d):
    raise RuntimeError("Database name not provided")
  return URL.create(drivername=DB_DRIVER, username=u, password=p,
                    database=d, host=DB_HOST, query={"unix_socket":s})

def connstr_from(args) -> str:
  """
  MySQL/MariaDB connection string based on the values of the args '<dbuser>',
  '<dbpass>', '<dbname>' and '<dbsocket>'
  """
  keys = ["<dbuser>", "<dbpass>", "<dbname>", "<dbsocket>"]
  return _connstr(*[args.get(k) for k in keys])

def connstr_env(varname) -> str:
  """
  Create a connection string from an env variable consisting of
  dbuser dbpass dbname and dbsocket, space-separated
  """
  return _connstr(*os.environ[varname].split(" "))
