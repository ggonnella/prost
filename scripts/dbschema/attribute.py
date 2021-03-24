#!/usr/bin/env python3
"""
DB Schema for tables storing attribute values.
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Text, Table, inspect, select
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.orm import declared_attr
from sqlalchemy_repr import PrettyRepresentableBase

Base = declarative_base(cls=PrettyRepresentableBase)

utf8_cs_args = {'mysql_charset': 'utf8', 'mysql_collate': 'utf8_bin'}

class AttributeDefinition(Base):
  """
  Describes an assembly attribute.

  If no definition is provided, the definition is given by the
  ontology link; otherwise the ontology link is a term which is
  related to the definition, in a way described in the definition.
  """
  __tablename__ = "pr_attribute_definition"
  id = Column(Integer, primary_key=True, autoincrement=False)
  name = Column(String(64), unique=True)
  datatype = Column(String(256), nullable=False)
  definition = Column(Text(1024))
  ontology_pfx = Column(String(64))
  ontology_term = Column(String(64))
  remark = Column(Text(1024))
  __table_args__ = utf8_cs_args

class AttributeValueMixin:
  accession = Column(String(64), primary_key = True)
  __table_args__ = utf8_cs_args
  __tablenamepfx__ = "pr_attribute_value_t"

  @declared_attr
  def __tablename__(cls):
    return cls.__tablenamepfx__ + str(cls.__tablenum__)

def _create_nth_table(n, engine):
  klass = type("AttributeValueT"+str(n), (AttributeValueMixin, Base),
               { "__tablenum__": n})
  klass.metadata.tables[klass.__tablename__].create(bind=engine.engine)

def nth_tablename(n):
  """
  Name of the n-th attribute_value table

  Note: the function does not check if a table with such name exists.
  """
  pfx = AttributeValueMixin.__tablenamepfx__
  return pfx+str(n)

def klass(n, engine):
  """
  Returns the class for the n-th table.
  """
  return Table(nth_tablename(n), Base.metadata, autoload_with=engine)

def high_table_num(engine):
  """
  Returns the highest table number of existing tables
  """
  insp = inspect(engine)
  tablenames = insp.get_table_names()
  pfx = AttributeValueMixin.__tablenamepfx__
  attrvaluetables = [tn for tn in tablenames if tn.startswith(pfx)]
  if len(attrvaluetables) == 0:
    return None
  tablenums = [int(tn[len(pfx):]) for tn in attrvaluetables]
  return max(tablenums)

def _attrnums_range_for_table(engine, tablename):
  cnames = [info["name"] for info in inspect(engine).get_columns(tablename)]
  attrnums = [int(cn.split("v")[0][len("attr"):]) for cn in \
        cnames if cn.startswith("attr") and not cn.endswith("c")]
  if attrnums:
    return min(attrnums), max(attrnums)
  else:
    return None, None

def nth_attrnums_range(engine, n):
  """
  Returns the range of value numbers of table with given table number
  """
  return _attrnums_range_for_table(engine, nth_tablename(n))

def _next_attrnum(engine, n):
  h = nth_attrnums_range(engine, n)[0]
  if h is None:
    if n == 0: return 0
    else: return _next_attrnum(engine, n-1)
  else: return h+1

TARGET_N_COLUMNS_FOR_TABLE = 64

def _n_columns_for_table(engine, tablename):
  return len(inspect(engine).get_columns(tablename))

def _new_attribute_coordinates(engine, ncols):
  """
  Return the table and attribute number for a new attribute
  for which ncols value columns and a computation id column are needed.
  Create the table if needed.

  A new table is created if the last created table has already value columns
  and the sum of the existing and new columns is higher than
  TARGET_N_COLUMNS_FOR_TABLE.
  """
  table_num = high_table_num(engine)
  attr_num = _next_attrnum(engine, table_num)
  nprevcols = _n_columns_for_table(engine, nth_tablename(table_num))
  if nprevcols > 1 and nprevcols + ncols + 1 > TARGET_N_COLUMNS_FOR_TABLE:
    table_num += 1
    _create_nth_table(table_num, engine)
  return table_num, attr_num

def _create_nth_attribute_columns(engine, datatypes, table_num, n):
  tablename = nth_tablename(table_num)
  if len(datatypes) == 1:
    coldefs = [f"attr{n}v {datatypes[0]}"]
  else:
    coldefs = [f"attr{n}v{i} {datatypes[i]}" for i in range(len(datatypes))]
  coldefs.append(f"attr{n}c {C_ID_TYPE}")
  coldefs = ", ".join(coldefs)
  engine.execute(f"ALTER TABLE {tablename} "+
                 f"ADD COLUMN ({coldefs})")

def table_num_for_nth_attr(engine, n):
  """
  Find number of table where the nth attribute is stored

  A linear search is performed. Null is returned if n is larger
  than the largest current attribute number.
  """
  assert(n >= 0)
  table_num = high_table_num(engine)
  l, h = nth_attrnums_range(engine, table_num)
  while h is None:
    if table_num == 0:
      return None
    else:
      table_num -= 1
      l, h = nth_attrnums_range(engine, table_num)
  if n > h:
    return None
  else:
    while n < l:
      table_num -= 1
      assert(table_num >= 0)
      l = nth_attrnums_range(engine, table_num)[0]
    return table_num

def _parse_datatype_def(datatype_def):
  result = []
  for elem in datatype_def.split(","):
    if elem.endswith("]"):
      dt, n = elem.split("[")
      n = int(n[:-1])
      assert(n > 1)
      for i in range(n):
        result.append(dt)
    else:
      result.append(elem)
  return result

def _get_table_numbers(engine):
  tablenames = inspect(engine).get_table_names()
  pfx = AttributeValueMixin.__tablenamepfx__
  attrvaluetables = [tn for tn in tablenames if tn.startswith(pfx)]
  return [int(tn[len(pfx):]) for tn in attrvaluetables]

def _get_attr_numbers(engine, n_tables):
  attrnums = []
  for tnum in range(n_tables):
    cnames = [info["name"] for info in \
        inspect(engine).get_columns(nth_tablename(tnum))]
    t_attrnums = set(int(cn.split("v")[0][len("attr"):]) for cn in \
        cnames if cn.startswith("attr") and not cn.endswith("c"))
    attrnums.append(t_attrnums)
  return attrnums

def _check_numbers(l, low=0):
  l.sort()
  prev = low-1
  for e in l:
    if e != prev+1:
      return False
    prev = e
  return True

C_ID_TYPE = "BINARY(16)"

def equiv(datatype1, datatype2):
  # not implemented, this should check if two datatype
  # definitions are equivalent, e.g. INT(11) == INTEGER
  return True

def _check_column(k, edt, tcols, checked_colnames, desc):
  if not k in tcols:
    raise ValueError(f"Missing column {k} ({desc})")
  checked_colnames.append(k)
  dt = tcols[k]
  if not equiv(dt, edt):
    raise ValueError(f"Wrong datatype for column {k} ({desc}): "+\
                     f"found {dt}, expected {edt}")

def _check_attribute_columns(engine, anum, tcols, datatypes, checked_colnames):
  k = f"attr{anum}c"
  _check_column(f"attr{anum}c", C_ID_TYPE, tcols, checked_colnames,
                f"computation ID of attribute {anum}")
  if len(datatypes) == 1:
    _check_column(f"attr{anum}v", datatypes[0], tcols, checked_colnames,
                  f"value column of attribute {anum}")
  else:
    for enum, edt in enumerate(datatypes):
      _check_column(f"attr{anum}v{enum}", edt, tcols, checked_colnames,
                    f"value element {enum} of attribute {anum}")

def check_consistency(engine):
  Session = sessionmaker(bind=engine)
  session = Session()
  tnums = _get_table_numbers(engine)
  if not _check_numbers(tnums):
    raise ValueError(f"The table numbers sequence is corrupted:\n{tnums}")
  anums = _get_attr_numbers(engine, len(tnums))
  expmin = 0
  allt_anums = []
  for tnum, t_anums in enumerate(anums):
    if t_anums:
      if not _check_numbers(list(t_anums), expmin):
        raise ValueError(f"The attribute numbers sequence is corrupted:\n{anums}")
      expmin = max(t_anums)+1
    allt_anums += list(t_anums)
  unexpected_adef = session.execute(select(AttributeDefinition).where(
                      AttributeDefinition.id.notin_(allt_anums))).all()
  if unexpected_adef:
    raise ValueError("Some attribute definitions do not correspond to "+
                     "columns in the attribute_value tables:\n"+
                     f"{unexpected_adef}")
  for tnum in tnums:
    tcols = {info[0]: info[1] for info in \
        engine.execute(f"DESCRIBE {nth_tablename(tnum)}").all()}
    checked_colnames = ["accession"]
    for anum in anums[tnum]:
      adef = session.execute(\
          select(AttributeDefinition).\
            where(AttributeDefinition.id==anum)).one()[0]
      datatypes = _parse_datatype_def(adef.datatype)
      _check_attribute_columns(engine, anum, tcols, datatypes, checked_colnames)
    unexpected_cols = set(tcols.keys())-set(checked_colnames)
    if unexpected_cols:
      raise ValueError(f"Unexpected columns in table {tnum}: {unexpected_cols}")

def create_attribute(engine, name, datatype_def, **kwargs):
  """
  Create a new attribute record in the attribute_definition table
  and reserve space in the attribute_values tables for storing
  the attribute values and computation IDs.

  Datatype definition, one of:
  - scalar: SQL datatype, including any parameter in ()
  - array: SQL datatype, followed by [<n>], with n integer > 0
  - struct: list of scalar and/or array, comma-sep, wo spaces

  e.g. TINYINT(1)[8],INTEGER,VARCHAR,CHAR(255)[2]
  """
  datatypes = _parse_datatype_def(datatype_def)
  tn, n = _new_attribute_coordinates(engine, len(datatypes))
  _create_nth_attribute_columns(engine, datatypes, tn, n)
  adef = AttributeDefinition(id = n, name = name,
                             datatype = datatype_def, **kwargs)
  Session = sessionmaker(bind=engine)
  session = Session()
  session.add(adef)
  session.commit()
