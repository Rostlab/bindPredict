""" Provides methods to parser output 
from secondary structure and solvent accessibility predictions like ProfPHD and Reprof.
"""

import sys
from Bio.Seq import Seq


class ProfParser(object):
    """ Parse and save information from ProfPDH or Reprof output.
    """

    # general
    seq = Seq("")

    # sec structure
    o_sec_struc = []
    p_sec_struc = []
    ri_sec_struc = []
    p_h = []
    p_e = []
    p_l = []
    o_t_h = []
    o_t_e = []
    o_t_l = []

    # solv access
    o_solv_acc = []
    p_solv_acc = []
    o_rel_solv_10 = []
    p_rel_solv_10 = []
    p_rel_solv = []
    p_10 = []
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

    # tmh
    omn = []
    pmn = []
    prmn = []
    ri_m = []
    p_m = []
    p_n = []
    o_t_m = []
    o_t_n = []

    def __init__(self, file_in):
        self.seq = Seq("")

        self.o_sec_struc = []
        self.p_sec_struc = []
        self.ri_sec_struc = []
        self.p_h = []
        self.p_e = []
        self.p_l = []
        self.o_t_h = []
        self.o_t_e = []
        self.o_t_l = []

        self.o_solv_acc = []
        self.p_solv_acc = []
        self.o_rel_solv_10 = []
        self.p_rel_solv_10 = []
        self.p_rel_solv = []
        self.p_10 = []
        self.ri_solv_acc = []
        self.o_be = []
        self.p_be = []
        self.o_bie = []
        self.p_bie = []
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

        self.omn = []
        self.pmn = []
        self.prmn = []
        self.ri_m = []
        self.p_m = []
        self.p_n = []
        self.o_t_m = []
        self.o_t_n = []

        self._parse_file(file_in)

    def _parse_file(self, file_in):
        rows = []

        with open(file_in) as f:
            for line in f:
                if '#' not in line:
                    row = line.strip().split("\t")
                    rows.append(row)

        tmp = rows.pop(0)

        col_map = {}
        for i in range(0, len(tmp)):
            col = tmp[i]
            col_map[col] = i

        self._parse_profphd(rows, col_map)

    def _parse_profphd(self, data, col_map):
        seq_str = ''

        for i in range(0, len(data)):
            el = data[i]
            aa = el[1]
            seq_str = seq_str + aa

            # get column indices
            o_sec_i = -1
            if 'OHEL' in col_map:
                o_sec_i = col_map['OHEL']
            p_sec_i = -1
            if 'PHEL' in col_map:
                p_sec_i = col_map['PHEL']
            ri_sec_i = -1
            if 'RI_S' in col_map:
                ri_sec_i = col_map['RI_S']
            p_h_i = -1
            if 'pH' in col_map:
                p_h_i = col_map['pH']
            p_e_i = -1
            if 'pE' in col_map:
                p_e_i = col_map['pE']
            p_l_i = -1
            if 'pL' in col_map:
                p_l_i = col_map['pL']
            o_t_h_i = -1
            if 'OtH' in col_map:
                o_t_h_i = col_map['OtH']
            o_t_e_i = -1
            if 'OtE' in col_map:
                o_t_e_i = col_map['OtE']
            o_t_l_i = -1
            if 'OtL' in col_map:
                o_t_l_i = col_map['OtL']

            o_solv_acc_i = -1
            if 'OACC' in col_map:
                o_solv_acc_i = col_map['OACC']
            p_solv_acc_i = -1
            if 'PACC' in col_map:
                p_solv_acc_i = col_map['PACC']
            o_rel_solv_10_i = -1
            if 'OREL' in col_map:
                o_rel_solv_10_i = col_map['OREL']
            p_rel_solv_10_i = -1
            if 'PREL' in col_map:
                p_rel_solv_10_i = col_map['PREL']
            p_10_i = -1
            if 'P10' in col_map:
                p_10_i = col_map['P10']
            ri_solv_i = -1
            if 'RI_A' in col_map:
                ri_solv_i = col_map['RI_A']
            o_be_i = -1
            if 'Obe' in col_map:
                o_be_i = col_map['Obe']
            p_be_i = -1
            if 'Pbe' in col_map:
                p_be_i = col_map['Pbe']
            o_bie_i = -1
            if 'Obie' in col_map:
                o_bie_i = col_map['Obie']
            p_bie_i = -1
            if 'Pbie' in col_map:
                p_bie_i = col_map['Pbie']
            o_t_0_i = -1
            if 'Ot0' in col_map:
                o_t_0_i = col_map['Ot0']
            o_t_1_i = -1
            if 'Ot1' in col_map:
                o_t_1_i = col_map['Ot1']
            o_t_2_i = -1
            if 'Ot2' in col_map:
                o_t_2_i = col_map['Ot2']
            o_t_3_i = -1
            if 'Ot3' in col_map:
                o_t_3_i = col_map['Ot3']
            o_t_4_i = -1
            if 'Ot4' in col_map:
                o_t_4_i = col_map['Ot4']
            o_t_5_i = -1
            if 'Ot5' in col_map:
                o_t_5_i = col_map['Ot5']
            o_t_6_i = -1
            if 'Ot6' in col_map:
                o_t_6_i = col_map['Ot6']
            o_t_7_i = -1
            if 'Ot7' in col_map:
                o_t_7_i = col_map['Ot7']
            o_t_8_i = -1
            if 'Ot8' in col_map:
                o_t_8_i = col_map['Ot8']
            o_t_9_i = -1
            if 'Ot9' in col_map:
                o_t_9_i = col_map['Ot9']

            omn_i = -1
            if 'OMN' in col_map:
                omn_i = col_map['OMN']
            pmn_i = -1
            if 'PMN' in col_map:
                pmn_i = col_map['PMN']
            prmn_i = -1
            if 'PRMN' in col_map:
                prmn_i = col_map['PRMN']
            ri_m_i = -1
            if 'RI_M' in col_map:
                ri_m_i = col_map['RI_M']
            p_m_i = -1
            if 'pM' in col_map:
                p_m_i = col_map['pM']
            p_n_i = -1
            if 'pN' in col_map:
                p_n_i = col_map['pN']
            o_t_m_i = -1
            if 'OtM' in col_map:
                o_t_m_i = col_map['OtM']
            o_t_n_i = -1
            if 'OtN' in col_map:
                o_t_n_i = col_map['OtN']

            if o_sec_i > -1:
                self.o_sec_struc.append(el[o_sec_i])
            if p_sec_i > -1:
                self.p_sec_struc.append(el[p_sec_i])
            if ri_sec_i > -1:
                self.ri_sec_struc.append(int(el[ri_sec_i]))
            if p_h_i > -1:
                self.p_h.append(int(el[p_h_i]))
            if p_e_i > -1:
                self.p_e.append(int(el[p_e_i]))
            if p_l_i > -1:
                self.p_l.append(int(el[p_l_i]))
            if o_t_h_i > -1:
                self.o_t_h.append(int(el[o_t_h_i]))
            if o_t_e_i > -1:
                self.o_t_e.append(int(el[o_t_e_i]))
            if o_t_l_i > -1:
                self.o_t_l.append(int(el[o_t_l_i]))

            if o_solv_acc_i > -1:
                self.o_solv_acc.append(int(el[o_solv_acc_i]))
            if p_solv_acc_i > -1:
                self.p_solv_acc.append(int(el[p_solv_acc_i]))
            if o_rel_solv_10_i > -1:
                self.o_rel_solv_10.append(int(el[o_rel_solv_10_i]))
            if p_rel_solv_10_i > -1:
                self.p_rel_solv_10.append(int(el[p_rel_solv_10_i]))
            if p_10_i > -1:
                self.p_10.append(int(el[p_10_i]))
            if ri_solv_i > -1:
                self.ri_solv_acc.append(int(el[ri_solv_i]))
            if o_be_i > -1:
                self.o_be.append(el[o_be_i])
            if p_be_i > -1:
                self.p_be.append(el[p_be_i])
            if o_bie_i > -1:
                self.o_bie.append(el[o_bie_i])
            if p_bie_i > -1:
                self.p_bie.append(el[p_bie_i])
            if o_t_0_i > -1:
                self.o_t_0.append(int(el[o_t_0_i]))
            if o_t_1_i > -1:
                self.o_t_1.append(int(el[o_t_1_i]))
            if o_t_2_i > -1:
                self.o_t_2.append(int(el[o_t_2_i]))
            if o_t_3_i > -1:
                self.o_t_3.append(int(el[o_t_3_i]))
            if o_t_4_i > -1:
                self.o_t_4.append(int(el[o_t_4_i]))
            if o_t_5_i > -1:
                self.o_t_5.append(int(el[o_t_5_i]))
            if o_t_6_i > -1:
                self.o_t_6.append(int(el[o_t_6_i]))
            if o_t_7_i > -1:
                self.o_t_7.append(int(el[o_t_7_i]))
            if o_t_8_i > -1:
                self.o_t_8.append(int(el[o_t_8_i]))
            if o_t_9_i > -1:
                self.o_t_9.append(int(el[o_t_9_i]))

            if omn_i > -1:
                self.omn.append(el[omn_i])
            if pmn_i > -1:
                self.pmn.append(el[pmn_i])
            if prmn_i > -1:
                self.prmn.append(el[prmn_i])
            if ri_m_i > -1:
                self.ri_m.append(el[ri_m_i])
            if p_m_i > -1:
                self.p_m.append(el[p_m_i])
            if p_n_i > -1:
                self.p_n.append(el[p_n_i])
            if o_t_m_i > -1:
                self.o_t_m.append(el[o_t_m_i])
            if o_t_n_i > -1:
                self.o_t_n.append(el[o_t_n_i])

        self.seq = Seq(seq_str)

    def _parse_reprof(self, data):
        seq_str = ''

        for i in range(0, len(data)):
            el = data[i]
            aa = el[1]
            seq_str = seq_str + aa
            self.p_sec_struc.append(el[2])
            self.ri_sec_struc.append(int(el[3]))
            self.p_h.append(int(el[4]))
            self.p_e.append(int(el[5]))
            self.p_l.append(int(el[6]))
            self.p_solv_acc.append(int(el[7]))
            self.p_rel_solv.append(int(el[8]))
            self.p_rel_solv_10.append(int(el[9]))
            self.ri_solv_acc.append(int(el[10]))
            self.p_be.append(el[11])
            self.p_bie.append(el[12])

        self.seq = Seq(seq_str)


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
