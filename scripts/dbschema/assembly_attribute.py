#!/usr/bin/env python3
"""
DB Schema for tables representing per-assembly attributes.
"""
from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Sequence, Column, Integer, Float, String, DateTime,\
                       Boolean, Text, Index
from sqlalchemy_repr import PrettyRepresentableBase
from sqlalchemy.dialects.mysql import DOUBLE

Base = declarative_base(cls=PrettyRepresentableBase)

utf8_cs_args = {'mysql_charset': 'utf8', 'mysql_collate': 'utf8_bin'}

class AssemblyIntValue(Base):
  __tablename__ = 'assembly_int_value'
  accession = Column(String(64), primary_key=True)
  attr_class = Column(String(64), primary_key=True, index=True)
  attr_instance = Column(String(64), primary_key=True)
  computation_id = Column(Integer)
  value = Column(Integer)
  __table_args__ = (Index('ix_assembly_int_value_attr',
                          'attr_class', 'attr_instance'), utf8_cs_args)

class AssemblyFloatValue(Base):
  __tablename__ = 'assembly_float_value'
  accession = Column(String(64), primary_key=True)
  attr_class = Column(String(64), primary_key=True, index=True)
  attr_instance = Column(String(64), primary_key=True)
  computation_id = Column(Integer)
  value = Column(DOUBLE)
  __table_args__ = (Index('ix_assembly_float_value_attr',
                          'attr_class', 'attr_instance'), utf8_cs_args)

