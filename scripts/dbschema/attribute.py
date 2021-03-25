#!/usr/bin/env python3
"""
DB Schema for tables storing attribute values.
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Text, Table,\
                       inspect, select, text
import sqlalchemy.types
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.orm import declared_attr
from sqlalchemy_repr import PrettyRepresentableBase
import ast
from collections import Counter

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
  name = Column(String(62), primary_key=True)
  datatype = Column(String(256), nullable=False)
  definition = Column(Text(4096))
  ontology_xref = Column(String(64))
  related_ontology_terms = Column(Text(4096))
  unit = Column(String(64))
  remark = Column(Text(4096))
  computation_group = Column(String(62), index=True)
  __table_args__ = utf8_cs_args

class AttributeValueMixin:
  accession = Column(String(64), primary_key = True)
  __table_args__ = utf8_cs_args
  __tablepfx__ = "pr_attribute_value_t"

  @declared_attr
  def __tablename__(cls):
    return cls.__tablepfx__ + str(cls.__tablesfx__)

class AttributeValueTables():

  @staticmethod
  def tablename(sfx):
    return AttributeValueMixin.__tablepfx__ + str(sfx)

  @staticmethod
  def normalize_suffix(sfx):
    return str(sfx).upper()

  @staticmethod
  def classname(sfx):
    return f"AttributeValueT{sfx}"

  def _init_attributes_maps(self):
    self._a2t = {}
    self._t2a = {}
    self._t2g = {}
    self._ncols = {}
    pfx = AttributeValueMixin.__tablepfx__
    inspector = inspect(self.engine)
    session = sessionmaker(bind=self.engine)()
    for tn in inspector.get_table_names():
      if tn.startswith(pfx):
        cnames = [c["name"] for c in inspector.get_columns(tn)]
        sfx = tn[len(pfx):]
        self._ncols[sfx] = len(cnames)
        cnames.remove("accession")
        anames = Counter(cn.rsplit("_", 1)[0] for cn in cnames \
            if not cn.endswith("_g"))
        self._t2a[sfx] = anames
        for an in anames:
          if an in self._a2t:
            raise RuntimeError(f"Attribute {an} found in multiple tables\n"+\
                               f"({self.tablename(self._a2t[an])} and {tn})")
          self._a2t[an] = sfx
        gnames = set(cn.rsplit("_", 1)[0] for cn in cnames \
            if cn.endswith("_g"))
        self._t2g[sfx] = {gname: set(session.execute(\
              select(AttributeDefinition.name).filter(\
            AttributeDefinition.computation_group == gname,
            AttributeDefinition.name.in_(anames))).scalars().all()) \
                for gname in gnames}

  def __init__(self, engine):
    self.engine = engine
    self._init_attributes_maps()

  @property
  def table_suffixes(self):
    return list(self._t2a.keys())

  @property
  def attribute_names(self):
    return list(self._a2t.keys())

  def table_for_attribute(self, attribute):
    return self.tablename(self._a2t[attribute])

  @staticmethod
  def tablesuffix(tablename):
    pfx = AttributeValueMixin.__tablepfx__
    if not tablename.startswith(pfx):
      raise RuntimeError(f"Table name does not start with prefix {pfx}")
    return tablename[len(pfx):]

  def table_attributes(self, tablename):
    return list(self._t2a[self.tablesuffix(tablename)].keys())

  def tables_for_attributes(self, attributes):
    result = {}
    for aname in attributes:
      if aname not in self._a2t:
        raise RuntimeError(f"Attribute not found: {aname}")
      tn = self.tablename(self._a2t[aname])
      if tn not in result: result[tn] = []
      result[tn].append(aname)
    return result

  def locations_for_attributes(self, attributes):
    """
    Computes the tablenames and column names where to store
    attribute values and computation IDs columns to set and to delete.

    If all elements of the group in a table are computed, the computation IDs
    to set is that of the group; in this case the computation IDs of the
    single elements must be set to NULL (to overwrite previous values, if any).

    If only part of the group is computed, the computation IDs to set
    are those of the attributes; the group ID is not set to NULL because
    it can still be valid for attributes which have not been recomputed.

    Return value:
      {"tables" =>
        {tablename => {"vcols": {attrname => vcols, ...},
                       "ccols_to_set": [...],
                       "ccols_to_unset": [...]},
          ...},
       "vcols" => []}
    """
    result = {"tables": {}, "vcols": []}
    for aname in attributes:
      if aname not in self._a2t:
        raise RuntimeError(f"Attribute not found: {aname}")
      t_sfx = self._a2t[aname]
      tn = self.tablename(t_sfx)
      if tn not in result["tables"]:
        result["tables"][tn] = {"vcols": {}, "ccols_to_set": [],
                                "ccols_to_unset": []}
      vcolnames = _vcolnames(aname, self._t2a[t_sfx][aname])
      result["tables"][tn]["vcols"][aname] = vcolnames
      result["vcols"] += vcolnames
    for tn, tdata in result:
      t_sfx = self.tablesuffix(tn)
      for g_name, g_members in self._t2g[t_sfx].items():
        if all(aname in result["tables"][tn]["vcols"] for aname in g_members):
          # whole group case
          result["tables"][tn]["ccols_to_set"].append(_gcolname(aname))
          result["tables"][tn]["ccols_to_unset"] += \
              [_ccolname(aname) for aname in g_members]
        else:
          result["tables"][tn]["ccols_to_set"] += \
              [_ccolname(aname) for aname in result["tables"][tn]["vcols"] \
                                if aname in g_members]
    return result

  def create_table(self, sfx):
    """
    Create new table with suffix <sfx>
    """
    sfx = self.normalize_suffix(sfx)
    if sfx in self._t2a:
      raise RuntimeError(f"Cannot create table: suffix {sfx} is not unique")
    klass = type(self.classname(sfx), (AttributeValueMixin, Base),
                 { "__tablesfx__": sfx})
    klass.metadata.tables[klass.__tablename__].create(bind=self.engine)
    self._t2a[sfx] = Counter()
    self._t2g[sfx] = {}
    self._ncols[sfx] = 1

  def get_class(self, sfx):
    """
    Class reflecting table with suffix <sfx>
    """
    return Table(self.tablename(self.normalize_suffix(sfx)), Base.metadata,
                 autoload_with=self.engine)

  def new_suffix(self) -> str:
    """
    Get a not-yet used numerical suffix
    """
    i = 0
    while True:
      if str(i) not in self._t2a:
        return str(i)
      i += 1

  #TARGET_N_COLUMNS = 64
  TARGET_N_COLUMNS = 6
  C_ID_TYPE = sqlalchemy.types.BINARY(16)

  def _place_for_new_attr(self, ncols, computation_group):
    """
    Table suffix for a new attribute for which ncols value columns and
    a computation id column are needed. Creates the table if needed.

    A new table is created if all tables have already value columns
    and the sum of the existing and new columns is higher than TARGET_N_COLUMNS.
    """
    for sfx in self.table_suffixes:
      if self._ncols[sfx] == 1:
        return sfx
      needed = ncols + 1
      if computation_group and computation_group not in self._t2g[sfx]:
        needed += 1
      if self._ncols[sfx] + needed <= self.TARGET_N_COLUMNS:
        return sfx
    sfx = self.new_suffix()
    self.create_table(sfx)
    return sfx

  @staticmethod
  def _vcolnames(a_name, nelems):
    if nelems == 1: return [f"{a_name}_v"]
    else:           return [f"{a_name}_v{i}" for i in range(nelems)]

  @staticmethod
  def _ccolname(a_name):
    return a_name+"_c"

  @staticmethod
  def _gcolname(grp_name):
    return grp_name+"_g"

  def _create_attr_cols(self, a_name, a_datatypes, g_name, t_sfx):
    if a_name in self._a2t:
      raise RuntimeError(f"Attribute {a_name} exists already in table"+\
                         self._a2t[a_name])
    colnames = self._vcolnames(a_name, len(a_datatypes))
    coldefs = [f"{cn} {dt}" for cn, dt in zip(colnames, a_datatypes)]
    coldefs.append(self._ccolname(a_name) + " " + str(self.C_ID_TYPE))
    if g_name and (g_name not in self._t2g[t_sfx]):
      coldefs.append(self._gcolname(g_name) + " " + str(self.C_ID_TYPE))
      self._t2g[t_sfx][g_name] = {a_name}
    coldefs_str = ", ".join(coldefs)
    tn = self.tablename(t_sfx)
    with self.engine.begin() as tr:
      tr.execute(text(f"ALTER TABLE {tn} ADD COLUMN ({coldefs_str})"))
    self._t2a[t_sfx][a_name] = len(a_datatypes)
    self._a2t[a_name] = t_sfx
    self._ncols[t_sfx] += len(coldefs)

  def _destroy_attr_cols(self, a_name, a_datatypes, g_name, t_sfx):
    colnames = self._vcolnames(a_name, len(a_datatypes))
    colnames.append(self._ccolname(a_name))
    tn = self.tablename(t_sfx)
    if g_name:
      self._t2g[t_sfx][g_name].remove(a_name)
      if len(self._t2g[t_sfx]) == 0:
        del self._t2g[t_sfx][g_name]
        colnames.append(self._gcolname(g_name))
    with self.engine.begin() as tr:
      dstr = ", ".join([f"DROP COLUMN {cn}" for cn in colnames])
      tr.execute(text(f"ALTER TABLE {tn} {dstr}"))
    del self._t2a[t_sfx][a_name]
    del self._a2t[a_name]
    self._ncols[t_sfx] -= len(colnames)

  @staticmethod
  def _parse_datatype_def(datatype_def):
    result = []
    for elem in datatype_def.split(";"):
      if elem.endswith("]"):
        dts, n = elem.split("[")
        n = int(n[:-1])
        assert(n > 1)
      else:
        dts = elem
        n = 1
      if dts.endswith(")"):
        dts_type, dts_sfx = elem.split("(")
        dts_params = [ast.literal_eval(p.strip()) \
                        for p in dts_sfx[:-1].split(",")]
      else:
        dts_type = dts
        dts_params = []
      if not hasattr(sqlalchemy.types, dts_type):
        raise ValueError("Unknown datatype in attribute datatype "+\
                         f"definition {dts_type}")
      dt = getattr(sqlalchemy.types, dts_type)(*dts_params)
      result += [dt]*n
    return result

  def create_attribute(self, name, datatype_def, **kwargs):
    """
    Create a new attribute record in the attribute_definition table
    and reserve space in the attribute_values tables for storing
    the attribute values and computation IDs.

    Datatype definition, one of:
    - scalar: datatype, including any parameter in ()
    - array: datatype, followed by [<n>], with n integer > 0
    - list of scalar and/or array, semicolon-sep, wo spaces

    Datatype is thereby any of the types defined
    in the sqlAlchemy.types module.

    e.g. Boolean[8];Integer;String(12);BINARY(16)[2]
    """
    a_datatypes = self._parse_datatype_def(datatype_def)
    t_sfx = self._place_for_new_attr(len(a_datatypes),
                                     kwargs.get("computation_group"))
    adef = AttributeDefinition(name = name,
                               datatype = datatype_def, **kwargs)
    session = sessionmaker(bind=self.engine)()
    session.add(adef)
    self._create_attr_cols(name, a_datatypes, adef.computation_group, t_sfx)
    session.commit()

  def destroy_attribute(self, name):
    session = sessionmaker(bind=self.engine)()
    adef=session.query(AttributeDefinition).filter(
        AttributeDefinition.name==name).first()
    if name not in self._a2t:
      raise RuntimeError(f"Attribute {name} does not exist")
    self._destroy_attr_cols(name, self._parse_datatype_def(adef.datatype),
                            adef.computation_group, self._a2t[name])
    session.delete(adef)
    session.commit()

  @staticmethod
  def _check_column(k, edt, cols, desc):
    if not k in cols:
      raise ValueError(f"Missing column {k} ({desc})")
    dt = cols[k]["type"]
    if str(dt) != str(edt):
      raise ValueError(f"Wrong datatype for column {k} ({desc}): "+\
                       f"found {dt}, expected {edt}")

  def check_consistency(self):
    session = sessionmaker(bind=self.engine)()
    inspector = inspect(self.engine)
    self._init_attributes_maps()
    unexpected_adef = session.execute(select(AttributeDefinition).where(
                      AttributeDefinition.name.notin_(
                        self.attribute_names))).all()
    if unexpected_adef:
      raise ValueError("Some attribute definitions do not correspond to "+
                       "columns in the attribute_value tables:\n"+
                       f"{unexpected_adef}")
    for sfx in self.table_suffixes:
      cols = {c["name"]: c for c in inspector.get_columns(self.tablename(sfx))}
      for aname in self._t2a[sfx]:
        adef = session.execute(\
          select(AttributeDefinition).\
            where(AttributeDefinition.name==aname)).scalars().one()
        self._check_column(self._ccolname(aname), self.C_ID_TYPE, cols,
                      f"computation ID of attribute {aname}")
        datatypes = self._parse_datatype_def(adef.datatype)
        if adef.computation_group:
          if adef.computation_group not in self._t2g[sfx]:
            raise RuntimeError("Column for computation group "+\
                f"{adef.computation_group} of attribute {aname} "+\
                f"not found in table {self.tablename(sfx)}")
        if len(datatypes) == 1:
          self._check_column(self._vcolnames(aname, 1)[0], datatypes[0], cols,
                        f"value column of attribute {aname}")
        else:
          vcolnames = self._vcolnames(aname, len(datatypes))
          for i, edt in enumerate(datatypes):
            self._check_column(vcolnames[i], edt, cols,
                          f"{i} element of value of attribute {aname}")
      for gname, group_members in self._t2g[sfx].items():
        if len(group_members) == 0:
          raise RuntimeError(f"Computation ID column for group {gname} "+\
              f"found in table {self.tablename(sfx)}, but no attribute "+\
              "of this group in the table.")
        self._check_column(self._gcolname(gname), self.C_ID_TYPE, cols,
                      f"computation ID of attribute group {gname}")
