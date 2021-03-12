from schema import Or, And, Use
import sys
import gzip

def open_maygz(fname):
  ft = open(fname, "rb")
  if ft.read(2) == b'\x1f\x8b':
    ft.close()
    return gzip.open(fname, "rt")
  else:
    ft.close()
    return open(fname)

open_file_or_stdin = Or(And(None, Use(lambda v: sys.stdin)), Use(open))

open_maygz_or_stdin = Or(And(None, Use(lambda v: sys.stdin)),
                         Use(open_maygz))

comments = Or(And(None, Use(lambda v: "#")), str)
delimiter = Or(And(None, Use(lambda v: "\t")), And(str, len))
