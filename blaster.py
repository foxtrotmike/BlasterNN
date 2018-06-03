#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 03:49:56 2018
An implementation of a BLAST-based nearest neighbor search
@author: fminhas
"""

from Bio import SeqIO
from Bio import SearchIO
import os
import numpy as np
from scipy.stats import rankdata
BLASTDir = "/mirror/pairpred_tools/ncbi-blast-2.2.30+-x64-linux/ncbi-blast-2.2.30+/bin/"
from joblib import Parallel, delayed

def blasterNN(testfile,trainfile,idx = "out"):
	"""
	Given a training and test FASTA file, this function will return a dictionary of the lowest e-values of a test protein against the training dataset. "idx" is used to assign an id to the output file.
	Note: it will creat a blast database file based on the training file so the location of the trainfile must be writeable. It will also generate the output file in the current directory.
	
	"""
	# Find the best matching protein (based on e-Value) in the training set using blast
    try:         
        cmd = BLASTDir+"makeblastdb -in "+trainfile+" -dbtype prot"
        os.system(cmd)
        cmd = BLASTDir+"blastp -db "+trainfile+" -query "+testfile+" -evalue 100 -out "+str(idx)+".pblast.txt"
        os.system(cmd)
        P = {}
        for qresult in SearchIO.parse(idx+".pblast.txt","blast-text"):
            if len(qresult):
                P[qresult.id]=np.min([hsp.evalue for hit in qresult for hsp in hit])
            else:
                P[qresult.id]=100.0
        return P
    except Exception:        
        print "Error Processing",idx
        return None
    
if __name__=='__main__':

	explist = [("testfile1.fasta,trfile1.fasta","1"),("testfile2.fasta,trfile2.fasta","2")]
	
    R = Parallel(n_jobs=4)(delayed(blasterNN)(testfile,trainfile,idx) for (test,train,idx) in explist)

