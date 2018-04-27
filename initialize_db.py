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
    type_strain = Column(Boolean)


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


class HmmHits(Base):
    __tablename__='HMMhits'
    hmmhitID = Column(Integer,primary_key=True)
    hmmhit = Column(Text,nullable=False) # row[3]
    orgname = Column(Text,ForeignKey('organisms.asmID')) #row[0].split(|)[0]
    seqid = Column(Integer,ForeignKey('Seqs.seqid')) #row[0].split(|)[1]
    hmmstart = Column(Integer) #row[15]
    hmmend = Column(Integer) #row[16]
    hmmlen = Column(Integer) #row[5]
    geneStart = Column(Integer) #row[19] (env)
    geneEnd = Column(Integer) #row[20] (env)
    genelen = Column(Integer) #row[2]
    evalue = Column(Float) #row[12] (i-Eval)
    score = Column(Float) #row[13]
    hmmcov = Column(Float)
    genecov = Column(Float)
    hmmLib = Column(Text)

class RnaHits(Base):
    __tablename__ = 'RNAHits'
    rnahit = Column(Text)
    orgname = Column(Text, ForeignKey('organisms.asmID'))
    seqid = Column(Integer,primary_key=True)
    eval = Column(Float)
    loc_start = Column(Integer)
    loc_end = Column(Integer)
    loc_strand = Column(Enum('+','-'))
    partialFlag = Column(Integer)
    seq = Column(Text, nullable=False)


engine = create_engine('sqlite:///uplb_mlst.db')
Base.metadata.drop_all(engine)
Base.metadata.create_all(engine)