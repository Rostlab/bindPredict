from config import FileSetter
from mmseqs_wrapper import MMSeqs2Wrapper

import math
from collections import defaultdict
import os


class MMseqsPredictor(object):

    def __init__(self, profile_db, blast_db, exclude_self_hits):
        self.profile_db = profile_db
        self.blast_db = blast_db
        self.exclude_self_hits = exclude_self_hits

    def get_mmseqs_hits(self, ids, evalue, criterion, set_name, hval=0):
        """
        Get MMseqs2 hits using local alignments at a specific E-value (H-value)
        :param ids: Ids to get hits for
        :param evalue: E-value cutoff to use for (1) MMseqs2 search and (2) filtering of hits
        :param criterion: eval or hval (should E-value or H-value be used to filter hits)
        :param set_name: Name for the query set; used for output files
        :param hval: H-value cutoff to use for filtering hits
        :return:
        """
        mmseqs_hits = defaultdict(defaultdict)

        set_ids = set(ids)
        # calc MMseqs2 alignments
        mmseqs_out = MMseqsPredictor.calc_mmseqs_results(evalue, set_name)

        if os.path.exists(mmseqs_out):
            # parse MMseqs2 output
            with open(mmseqs_out) as read_mmseqs:
                for line in read_mmseqs:
                    splitted_line = line.strip().split()
                    i = splitted_line[0]
                    if i in set_ids:
                        hit = splitted_line[1]
                        if i != hit or not self.exclude_self_hits:
                            e_val = float(splitted_line[2])

                            nident = int(splitted_line[3])  # number of identical positions
                            nmismatch = int(splitted_line[4])  # number of mismatches

                            lali = nident + nmismatch
                            pide = nident / lali

                            h_val = MMseqsPredictor._calc_hssp(pide, lali)

                            if criterion == 'eval' and e_val <= float(evalue):
                                # filter hits by E-value
                                mmseqs_hits[i][hit] = {'pide': pide, 'lali': lali, 'eval': e_val, 'hval': h_val}
                            elif criterion == 'hval' and h_val >= float(hval):
                                # filter hits by H-value
                                mmseqs_hits[i][hit] = {'pide': pide, 'lali': lali, 'eval': e_val, 'hval': h_val}

        # determine hit with minimal E-value/maximum H-value
        best_hits = defaultdict(dict)
        for i in mmseqs_hits.keys():
            max_hit = ''
            min_eval = 1e5
            max_pide = 0
            max_hval = -100

            for h in mmseqs_hits[i].keys():
                e_val = mmseqs_hits[i][h]['eval']
                pide = mmseqs_hits[i][h]['pide']
                h_val = mmseqs_hits[i][h]['hval']

                if criterion == 'eval' and (e_val < min_eval or (e_val <= min_eval and pide > max_pide)):
                    max_hit = h
                    min_eval = e_val
                    max_pide = pide
                if criterion == 'hval' and (h_val > max_hval or (h_val >= max_hval and pide > max_pide)):
                    max_hit = h
                    max_hval = h_val
                    max_pide = pide

            if criterion == 'eval':
                best_hits[i][max_hit] = min_eval
            else:
                best_hits[i][max_hit] = max_hval

        return best_hits

    @staticmethod
    def calc_mmseqs_results(evalue, set_name):
        """
        Calculate MMseqs2 alignments using profiles against big80
        :param evalue: e-value threshold for final alignments
        :param set_name: Meaningful set name to ensure that previous results are not overwritten/re-used
        :return: File to which MMseqs2 results were written
        """

        fasta_in = FileSetter.query_set()
        aln_out = '{}/mmseqs_big80_results_{}.m8'.format(FileSetter.mmseqs_output(), set_name)

        # mmseqs profile search against big80
        profile_db = FileSetter.profile_db()  # use pre-computed database from big80
        lookup_db = FileSetter.lookup_db()  # use pre-computed database of lookup set

        mmseqs = MMSeqs2Wrapper(FileSetter.mmseqs_path(), FileSetter.mmseqs_output(), set_name)

        # create database of the query set
        query_db = mmseqs.create_db(fasta_in, set_name)

        # create profiles against big80
        aln_res = mmseqs.seq2seq_search(query_db, profile_db, num_iterations=2)
        query_profiles = mmseqs.alnres2prof(query_db, profile_db, aln_res)

        # search profiles against lookup db
        profile_hits_raw = mmseqs.prof2seq_search(query_profiles, lookup_db, 0, evalue)
        mmseqs.convertali(query_db, lookup_db, profile_hits_raw, aln_out)

        return aln_out

    @staticmethod
    def _calc_hssp(pide, lali):
        """
        Calculate H-value
        :param pide: Percentage pairwise sequence identity
        :param lali: Alignment length (excl. gaps)
        :return:
        """

        pide = pide * 100

        if lali <= 11:
            return pide - 100
        elif lali > 450:
            return pide - 19.5
        else:
            exp = -0.32 * (1 + math.exp(-lali / 1000))
            return pide - (480 * (lali ** exp))
