#!/usr/bin/env python3
"""
DB Schema for tables storing attribute values.
"""
from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Sequence, Column, Integer, Float, String, DateTime,\
                       Boolean, Text, Index
from sqlalchemy_repr import PrettyRepresentableBase
from sqlalchemy.dialects.mysql import DOUBLE, BINARY
from sqlalchemy.orm import declared_attr

Base = declarative_base(cls=PrettyRepresentableBase)

utf8_cs_args = {'mysql_charset': 'utf8', 'mysql_collate': 'utf8_bin'}

def tablename(valuetype, scored):
  parts = ["pr", valuetype]
  if scored: parts.append("scored")
  parts.append("attribute_value")
  return "_".join(parts)

# defines an order of the columns in data files resulting from
# a batch computation (therefore computation_id is left out)
def datafile_columns(scored):
  result = ["accession", "attr_class", "attr_instance", "value"]
  if scored: result.append("score")
  return result

__valuecoltype__ = {"int": Integer, "float": DOUBLE}

class AttributeValueMixin:
  accession = Column(String(64), primary_key = True)
  attr_class = Column(String(64), primary_key = True, index=True)
  attr_instance = Column(String(64), primary_key = True)
  computation_id = Column(BINARY(16))

  @declared_attr
  def value(cls):
    return Column(__valuecoltype__[cls.__valuetype__])

  @declared_attr
  def __tablename__(cls):
    return tablename(cls.__valuetype__, hasattr(cls, "score"))

  @declared_attr
  def __table_args__(cls):
    return (Index("_".join(["ix", cls.__tablename__, "attr"]),
                    "attr_class", "attr_instance"), utf8_cs_args)

class ScoredAttributeValueMixin(AttributeValueMixin):
  score = Column(DOUBLE)

class IntAttributeValue(AttributeValueMixin, Base):
  __valuetype__ = 'int'

class FloatAttributeValue(AttributeValueMixin, Base):
  __valuetype__ = 'float'

class IntScoredAttributeValue(ScoredAttributeValueMixin, Base):
  __valuetype__ = 'int'

class FloatScoredAttributeValue(ScoredAttributeValueMixin, Base):
  __valuetype__ = 'float'

