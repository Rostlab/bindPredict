""" Provides methods to parser output 
from secondary structure and solvent accessibility predictions like ProfPHD and Reprof.
"""

from Bio.Seq import Seq


class ProfParser:
    """ Parse and save information from ProfPDH or Reprof output.
    """

    seq = Seq("")
    o_sec_struc = []
    p_sec_struc = []
    ri_sec_struc = []
    p_h = []
    p_e = []
    p_l = []
    o_t_h = []
    o_t_e = []
    o_t_l = []
    o_solv_acc = []
    p_solv_acc = []
    o_rel_solv_10 = []
    p_rel_solv_10 = []
    p_rel_solv = []
    ri_solv_acc = []
    o_be = []
    p_be = []
    o_bie = []
    p_bie = []
    o_t_0 = []
    o_t_1 = []
    o_t_2 = []
    o_t_3 = []
    o_t_4 = []
    o_t_5 = []
    o_t_6 = []
    o_t_7 = []
    o_t_8 = []
    o_t_9 = []

    def __init__(self, file_in, method):
        self._parse_file(file_in, method)

    def _parse_file(self, file_in, method):
        rows = []

        with open(file_in) as f:
            for line in f:
                if '#' not in line:
                    row = line.strip().split("\t")
                    rows.append(row)

        rows.pop(0)

        if method == 'reprof':
            self._parse_reprof(rows)
        elif method == 'profphd':
            self._parse_profphd(rows)
        else:
            raise AssertionError('No valid method specified')

    def _parse_profphd(self, data):
        seq_str = ''
        self.p_rel_solv = []

        for i in range(0, len(data) - 1):
            el = data[i]
            aa = el[1]
            seq_str = seq_str + aa
            self.o_sec_struc.append(el[2])
            self.p_sec_struc.append(el[3])
            self.ri_sec_struc.append(el[4])
            self.o_solv_acc.append(el[5])
            self.p_solv_acc.append(el[6])
            self.o_rel_solv_10.append(el[7])
            self.p_rel_solv_10.append(el[8])
            self.ri_solv_acc.append(el[9])
            self.p_h.append(el[10])
            self.p_e.append(el[11])
            self.p_l.append(el[12])
            self.o_be.append(el[13])
            self.p_be.append(el[14])
            self.o_bie.append(el[15])
            self.p_bie.append(el[16])
            self.o_t_h.append(el[17])
            self.o_t_e.append(el[18])
            self.o_t_l.append(el[19])
            self.o_t_0.append(el[20])
            self.o_t_1.append(el[21])
            self.o_t_2.append(el[22])
            self.o_t_3.append(el[23])
            self.o_t_4.append(el[24])
            self.o_t_5.append(el[25])
            self.o_t_6.append(el[26])
            self.o_t_7.append(el[27])
            self.o_t_8.append(el[28])
            self.o_t_9.append(el[29])

        self.seq = Seq(seq_str)

    def _parse_reprof(self, data):
        seq_str = ''

        self.o_sec_struc = []
        self.o_solv_acc = []
        self.o_rel_solv_10 = []
        self.o_be = []
        self.o_bie = []
        self.o_t_h = []
        self.o_t_e = []
        self.o_t_l = []
        self.o_t_0 = []
        self.o_t_1 = []
        self.o_t_2 = []
        self.o_t_3 = []
        self.o_t_4 = []
        self.o_t_5 = []
        self.o_t_6 = []
        self.o_t_7 = []
        self.o_t_8 = []
        self.o_t_9 = []

        for i in range(0, len(data) - 1):
            el = data[i]
            aa = el[1]
            seq_str = seq_str + aa
            self.p_sec_struc.append(el[2])
            self.ri_sec_struc.append(el[3])
            self.p_h.append(el[4])
            self.p_e.append(el[5])
            self.p_l.append(el[6])
            self.p_solv_acc.append(el[7])
            self.p_rel_solv.append(el[8])
            self.p_rel_solv_10.append(el[9])
            self.ri_solv_acc.append(el[10])
            self.p_be.append(el[11])
            self.p_bie.append(el[12])

        self.seq = Seq(seq_str)
