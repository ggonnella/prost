from contextlib import contextmanager
import MySQLdb

def connect(default_args):
  return MySQLdb.connect(host="localhost",
                         user=args["<dbuser>"],
                         passwd=args["<dbpass>"],
                         db=args["<dbname>"],
                         unix_socket=args["<dbsocket>"],
                         use_unicode=True)

@contextmanager
def cursor_from(default_args):
  db = connect(default_args)
  cursor = db.cursor()
  try:
    yield cursor
  finally:
    cursor.close()
  db.commit()

def connect_and_execute(default_args, statements):
  with cursor_from(default_args) as c:
    for statement in statements:
      c.execute(statement)
