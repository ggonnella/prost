from schema import Or, And, Use
import sys

open_file_or_stdin = Or(And(None, Use(lambda v: sys.stdin)), Use(open))
comments = Or(And(None, Use(lambda v: "#")), str)
delimiter = Or(And(None, Use(lambda v: "\t")), And(str, len))
