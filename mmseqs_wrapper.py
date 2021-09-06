#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 11:34:00 2021

@author: mheinzinger
"""

import subprocess
import os

"""
Overview over pipeline steps:
    1. Create DBs for queries and lookup
        mmseqs createdb test_set.fasta test_set
        mmseqs createdb lookup_set.fasta lookup_seq
    2. Use sequence-sequence search to search test against lookup
        mmseqs search test_set lookup_set alnResults tmp --num-iterations 3 -a
    3. Compute profiles from hits
        mmseqs result2profile test_set lookup_set alnResults test_set_profile 
    4. Search test set profiles against lookup sequences to detect remote hits
        mmseqs search test_set_profile lookup_set alnResults tmp --min-seq-id 0.2 -s 7.5 --max-seqs 10000
"""


class MMSeqs2Wrapper(object):
    
    def __init__(self, mmseqs_path, work_dir, query_set):
        self.mmseqs = mmseqs_path
        self.work_dir = work_dir  # Directory for storing final results
        self.tmp_dir = '{}/tmp'.format(self.work_dir)  # Directory for storing intermediate results
        self.query_set = query_set
        
    def create_db(self, fasta_p, set_name):
        db_res = '{}/{}'.format(self.tmp_dir, set_name)
        if not os.path.exists(db_res):
            subprocess.call([str(self.mmseqs), "createdb", fasta_p, db_res])
        return db_res
    
    def seq2seq_search(self, query_db, lookup_db, num_iterations=3):
        aln_res = '{}/seq2seq_alnRes_{}'.format(self.tmp_dir, self.query_set)

        subprocess.call([str(self.mmseqs), "search", query_db, lookup_db, aln_res, str(self.tmp_dir),
                        "--num-iterations", str(num_iterations), "-a"])
        return aln_res
    
    def alnres2prof(self, query_db, lookup_db, aln_res):
        query_profiles = '{}/test_set_profiles_{}'.format(self.tmp_dir, self.query_set)

        subprocess.call([str(self.mmseqs), "result2profile", query_db, lookup_db, aln_res, query_profiles])
        return query_profiles
    
    def prof2seq_search(self, query_profiles, lookup_db, min_seq_id, e_val, sensitivity=7.5, max_seqs=100000):
        aln_res = '{}/prof2seq_alnRes_{}'.format(self.tmp_dir, self.query_set)

        subprocess.call([str(self.mmseqs), "search", query_profiles, lookup_db, aln_res, str(self.tmp_dir),
                         "--min-seq-id", str(min_seq_id), "-s", str(sensitivity), "--max-seqs", str(max_seqs),
                         "-e", e_val, "-a"])
        return aln_res

    def convertali(self, query_db, lookup_db, profile_hits_raw, output_file):
        subprocess.call([str(self.mmseqs), "convertalis", query_db, lookup_db, profile_hits_raw, output_file,
                         "--format-output",  "query,target,evalue,nident,mismatch,qstart,tstart,qaln,taln"])
        return output_file
