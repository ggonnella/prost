#!/usr/bin/env python3
"""
Input: TSV containing at least two columns:
  (a) JSON information
  (b) an ID / primary key

Output: TSV containing a line for each leaf of the JSON tree
        content of the lines:
  (a) path to the leaf
  (b) leaf key
  (c) record primary key
  (d) leaf value (as JSON)

Explanation: JSON is seen by the script as a tree with
named nodes. Names are not unique. If a node contains:
  - a dict, then it is an internal tree node:
    - child nodes are the dict values
    - their names are the keys
  - a list, then it is an internal tree node:
    - child nodes are the list elements
    - their names are <list_name>[]
  - a scalar, then it is a leaf:
    - node name is either "/" (root), a dictionary key
      or list element index (see above)
    - the scalar is the leaf content

Usage:
  json_tree_leaves.py [options] <tsvfile> <jsoncol> <idcol>

Arguments:
  tsvfile:  the input TSV file
  jsoncol:  1-based column num, column, JSON information
  idcol:    1-based column num, IDs / primary keys

Options:
{common}
"""

from docopt import docopt
from schema import And, Use
from lib import scripts, snake
import json

def print_leaves(subtree, primarykey, pfx, last):
  if isinstance(subtree, dict):
    for k, v in subtree.items():
      print_leaves(v, primarykey, pfx + "\\" + k, k)
  elif isinstance(subtree, list):
    for v in subtree:
      print_leaves(v, primarykey, pfx+"[]", last+"[]")
  else:
    print("\t".join([pfx, last, primarykey,
           json.dumps(subtree)]))

def main(arguments):
  for line in arguments["<tsvfile>"]:
    elems = line.rstrip().split("\t")
    primarykey = elems[arguments["<idcol>"]]
    tree = json.loads(elems[arguments["<jsoncol>"]])
    print_leaves(tree, primarykey, "", "")

def validated(args):
  colnum = And(Use(lambda i: int(i)-1), lambda i: i>=0)
  return scripts.validate(args, {"<tsvfile>": Use(open),
           "<jsoncol>": colnum, "<idcol>": colnum})

if "snakemake" in globals():
  args = snake.args(snakemake,
      input=["<tsvfile>"], params=["<jsoncol>", "<idcol>"])
  main(validated(args))
elif __name__ == "__main__":
  args = docopt(__doc__.format(common=scripts.args_doc), version=0.1)
  main(validated(args))
