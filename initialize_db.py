#import os,sys
from sqlalchemy import Column,ForeignKey,Integer,String,Text,Enum,Float,Boolean
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine

Base = declarative_base()

class Organism(Base):
    __tablename__ = 'organisms'
    id = Column(Integer,primary_key=True,nullable=False)
    asmID = Column(Text,nullable=False,unique=True)
    speciesName = Column(Text)

    taxID = Column(Integer)
    ncbiID = Column(Integer)
    rsCat = Column(Text)
    anomaly= Column(Text)
    asmStatus = Column(Text)
    asmVersion = Column(Integer,nullable=False)
    type_strain = relationship("Phylogeny",back_populates='type_strain')

    cds = relationship("CDS",back_populates='seqid')
    hmmHits = relationship("HMMhits",back_populates='hmmhit')
    rnaHits = relationship("RNAHits",back_populates='rnahit')


class Phylogeny(Base):
    __tablename__ = 'phylogeny'
    asmID = Column(Text,ForeignKey('organisms.asmID'))
    familyID = Column(Integer)
    familyName = Column(Text)
    genusName = Column(Text)
    genusID = Column(Integer)
    orderName = Column(Text)
    orderID = Column(Integer)
    phylumID = Column(Integer)
    phylumName = Column(Text)
    speciesName = Column(Text)
    type_strain = Column(Boolean)
    taxID = Column(Integer,primary_key=True)



class CDS(Base):
    __tablename__= 'Seqs'
    seqid = Column(Integer,primary_key=True,nullable=False)
    # naming inconsistency to maintain autoMLST compatibility
    orgname = Column(Text,ForeignKey('organisms.asmID'))
    gene = Column(Text)
    description = Column(Text)
    source = Column(Text)
    loc_start = Column(Integer)
    loc_end = Column(Integer)
    loc_strand = Column(Enum('+','-'))
    acc = Column(Text)
    lastscan = Column(Integer)
    naseq = Column(Text)
    aaseq = Column(Text,nullable=False)

    hmmHits = relationship("HMMhits",back_populates='hmmhit')

class HmmHits(Base):
    __tablename__='HMMhits'
    hmmhitID = Column(Integer,primary_key=True)
    hmmhit = Column(Text,nullable=False)
    orgname = Column(Text,ForeignKey('organisms.name'))
    seqid = Column(Integer,ForeignKey('Seqs.seqid'))
    hmmstart = Column(Integer)
    hmmend = Column(Integer)
    hmmlen = Column(Integer)
    geneStart = Column(Integer)
    geneEnd = Column(Integer)
    genelen = Column(Integer)
    evalue = Column(Float)
    score = Column(Float)
    bias = Column(Float)
    iscore = Column(Float)
    flags = Column(Integer)
    hmmcov = Column(Float)
    genecov = Column(Float)
    hmmLib = Column(Text)

    cds = relationship("Seqs",back_populates='orgname')

class RnaHits(Base):
    __tablename__ = 'RNAHits'
    rnahit = Column(Text)
    orgname = Column(Text, ForeignKey('organisms.name'))
    seqid = Column(Integer,primary_key=True)
    score = Column(Float)
    loc_start = Column(Integer)
    loc_end = Column(Integer)
    loc_strand = Column(Enum('+','-'))
    seq = Column(Text, nullable=False)

if __name__ == '__main__':
    engine = create_engine('sqlite:///mlst_sql.db')
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)