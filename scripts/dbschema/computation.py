#!/usr/bin/env python3
"""
DB Schema for tables representing computations.

Tables storing the results of computations can refer to this table for the sake
of data traceability and computation reproducibility.
"""
from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Sequence, Column, Integer, Float, String, DateTime,\
                       Boolean, Text, Enum
from sqlalchemy_repr import PrettyRepresentableBase
from sqlalchemy.dialects.mysql import DOUBLE

Base = declarative_base(cls=PrettyRepresentableBase)

utf8_cs_args = {'mysql_charset': 'utf8', 'mysql_collate': 'utf8_bin'}

class Computation(Base):
  __tablename__ = 'computation'
  id = Column(Integer, primary_key=True, autoincrement=True)
  plugin = Column(String(256), nullable=False)
  version = Column(String(64), nullable=False)
  parameters = Column(Text(4096))
  unit = Column(Enum("assembly"), nullable=False, default="assembly")
  n_units = Column(Integer, nullable=False)
  reason = Column(Enum("new_data", "new_attributes", "improve_precision"),
                  nullable=False)
  comp_status = Column(Enum("running", "completed", "partial", "aborted"),
                  nullable=False, default="running")
  system_id = Column(String(64))
  user_id = Column(String(64))
  time_start = Column(DateTime)
  time_end = Column(DateTime)
  used_resources = Column(Text(4096))
  __table_args__ = utf8_cs_args

