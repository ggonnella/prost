import csv
import sys

def decomment(f, pfx, verbose=False):
  for line in f:
    if not pfx or not line.startswith(pfx):
      yield line
    elif verbose:
      sys.stderr.write("# skipped internal comment line: "+line)

def get_fieldnames_for(fieldslist, delimiter, commentspfx, verbose, f):
  if fieldslist:
    fieldnames=fieldslist.split(",")
  else:
    line = f.readline()
    while delimiter not in line:
      if verbose:
        sys.stderr.write("# skipped initial comment line: "+line)
      line = f.readline()
    if line.startswith(commentspfx):
      line = line[len(commentspfx):]
    fieldnames = [x.strip() for x in line.split(delimiter)]
  if verbose:
    sys.stderr.write("# table field names: "+", ".join(fieldnames)+"\n")
  return fieldnames

def get_dict_reader_for(fieldnames, delimiter, commentspfx, verbose, f):
  return csv.DictReader(
      decomment(f, commentspfx, verbose),
      delimiter=delimiter,
      fieldnames=fieldnames,
      quoting=csv.QUOTE_NONE)

def get_fieldnames(args, f):
  return get_fieldnames_for(args["--fields"], args["--delimiter"],
                        args["--comments"], args["--verbose"], f)

def get_dict_reader(args, f):
  return get_dict_reader_for(get_fieldnames(args, f),
      args["--delimiter"], args["--comments"], args["--verbose"], f)

def n_header_lines(filename: str, pfx: str = "#") -> int:
  """Number of header lines, identified by a given prefix

  Args:
    filename: data of the table file
    pfx:      header lines prefix, defaults to "#"

  Returns:
    0 if pfx is empty, otherwise number of initial lines starting with pfx
  """
  result = 0
  if not pfx:
    return 0
  with open(filename) as f:
    for line in f:
      if line.startswith(pfx):
        result += 1
      else:
        return result


