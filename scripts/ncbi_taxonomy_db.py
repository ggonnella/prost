#!/usr/bin/env python3

from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Sequence, Column, Integer, String, \
                       Index, DateTime, Text, Boolean
from sqlalchemy.orm import relationship

# Nt = NCBI Taxonomy

Base = declarative_base()

utf8_cs_args = {'mysql_charset': 'utf8', 'mysql_collate': 'utf8_bin'}

class NtName(Base):
  __tablename__ = 'nt_names'
  id = Column(Integer, Sequence('nt_name_id_seq'), primary_key=True)
  time_created = Column(DateTime, server_default=func.now())
  tax_id = Column(Integer)
  name_txt = Column(String(256))
  unique_name = Column(String(256))
  name_class = Column(String(64))
  __table_args__ = utf8_cs_args

  @classmethod
  def file_column_names(cls):
    return ["tax_id", "name_txt", "unique_name", "name_class"]

  @classmethod
  def create_indices(cls, engine):
    indices = [Index('nt_tax_id_index', cls.tax_id),
               Index('nt_names_txt_index', cls.name_txt),
               Index('nt_names_row_index', cls.tax_id, cls.name_txt,
                     cls.unique_name, cls.name_class, unique=True)]
    for i in indices:
      i.create(bind=engine)

class NtGencode(Base):
  __tablename__ = 'nt_gencode'
  genetic_code_id = Column(Integer, primary_key=True, autoincrement=False)
  time_created = Column(DateTime, server_default=func.now())
  abbreviation = Column(String(64))
  name = Column(String(256))
  cde = Column(String(64))
  starts = Column(String(64))
  __table_args__ = utf8_cs_args

  @classmethod
  def file_column_names(cls):
    return ["genetic_code_id", "abbreviation", "name", "cde", "starts"]

  @classmethod
  def create_indices(cls, engine):
    pass

class NtMerged(Base):
  __tablename__ = 'nt_merged'
  id = Column(Integer, Sequence('nt_merged_id_seq'), primary_key=True)
  time_created = Column(DateTime, server_default=func.now())
  old_tax_id = Column(Integer)
  new_tax_id = Column(Integer)
  __table_args__ = utf8_cs_args

  @classmethod
  def file_column_names(cls):
    return ["old_tax_id", "new_tax_id"]

  @classmethod
  def create_indices(cls, engine):
    indices = [Index('nt_merged_old_tax_id_index', cls.old_tax_id),
               Index('nt_merged_new_tax_id_index', cls.new_tax_id)]
    for i in indices:
      i.create(bind=engine)

class NtDivision(Base):
  __tablename__ = 'nt_division'
  division_id = Column(Integer, primary_key=True, autoincrement=False)
  time_created = Column(DateTime, server_default=func.now())
  division_cde = Column(String(3))
  division_name = Column(String(64))
  comments = Column(String(256))
  __table_args__ = utf8_cs_args

  @classmethod
  def file_column_names(cls):
    return ["division_id", "division_cde", "division_name", "comments"]

  @classmethod
  def create_indices(cls, engine):
    pass

class NtDelnode(Base):
  __tablename__ = 'nt_delnodes'
  tax_id = Column(Integer, primary_key=True, autoincrement=False)
  __table_args__ = utf8_cs_args

  @classmethod
  def file_column_names(cls):
    return ["tax_id"]

  @classmethod
  def create_indices(cls, engine):
    pass

class NtCitation(Base):
  __tablename__ = 'nt_citations'
  cit_id = Column(Integer, primary_key=True, autoincrement=False)
  cit_key = Column(String(256))
  pubmed_id = Column(Integer)
  medline_id = Column(Integer)
  url = Column(String(1024))
  text = Column(Text)
  taxid_list = Column(Text(1000000))
  __table_args__ = utf8_cs_args

  @classmethod
  def file_column_names(cls):
    return ["cit_id", "cit_key", "pubmed_id", "medline_id", "url", "text",
            "taxid_list"]

  @classmethod
  def create_indices(cls, engine):
    pass

class NtNode(Base):
  __tablename__ = 'nt_nodes'
  tax_id = Column(Integer, primary_key=True, autoincrement=False)
  parent_tax_id = Column(Integer)
  rank = Column(String(64))
  embl_code = Column(String(64))
  division_id = Column(Integer)
  inherited_div_flag = Column(Boolean)
  genetic_code_id = Column(Integer)
  inherited_GC_flag = Column(Boolean)
  mitochondrial_genetic_code_id = Column(Integer)
  inherited_MGC_flag = Column(Boolean)
  GenBank_hidden_flag = Column(Boolean)
  hidden_subtree_root_flat = Column(Boolean)
  comments = Column(String(1024))
  __table_args__ = utf8_cs_args

  @classmethod
  def file_column_names(cls):
    return ["tax_id", "parent_tax_id", "rank", "embl_code", "division_id",
            "inherited_div_flag", "genetic_code_id", "inherited_GC_flag",
            "mitochondrial_genetic_code_id", "inherited_MGC_flag",
            "GenBank_hidden_flag", "hidden_subtree_root_flat", "comments"]

  @classmethod
  def create_indices(cls, engine):
    indices = [Index('nt_nodes_parent_tax_id_index', cls.parent_tax_id),
               Index('nt_nodes_rank_index', cls.rank)]
    for i in indices:
      i.create(bind=engine)

tablename2class = {
  'nt_names': NtName,
  'nt_gencode': NtGencode,
  'nt_merged': NtMerged,
  'nt_division': NtDivision,
  'nt_delnodes': NtDelnode,
  'nt_citations': NtCitation,
  'nt_nodes': NtNode,
}

