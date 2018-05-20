from Bio.Seq import Seq
from collections import defaultdict
import re
import os.path
import numpy
from scripts.parser import ProfParser

class Protein(object):

    """ A protein in the context of binding site prediction """

    def __init__(self):
        self.snap_file = None
        self.evc_file = None
        self.dist_file = None
        self.evmut_file = None
        self.cons_file = None
        self.solv_file = None
        self.length = 0
        self.length_cov = 0
        self.predictions = dict()
        self.performance = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0, 'Cov': 0, 'Prec': 0, 'F1':0, 'Acc':0}
        self.scores = defaultdict(list)
        self.binding_res = []
        self.seq_dist = 5
        self.threshold = 0.0
        self.average = 0.0
        self.factor = 2
        self.cs_dist_thresh = 8.0
        self.cc_dist_thresh = 20.0
        self.blosum_file = None
        self.blosum_thresh = 0.0
        self.blosum_mat = defaultdict(dict)
        self.snap_percentage = 0.3
        self.evc_percentage = 0.4
        self.abs_vals = {'A':106, 'C':135, 'F':197, 'I':169, 'M':188, 'Q':198, 'T':142, 'X':180, 'Y':222,
			 'B':160, 'D':163, 'G':84, 'K':205, 'N':157, 'R':248, 'V':142, 'Z':196, 'E':194,
			 'H':184, 'L':164, 'P':136, 'S':130, 'W':227, 'U':135}


    def calc_scores(self):
        """ Calculate cum scores, clustering coefficients, SNAP2/BLOSUM, EVmut/BLOSUM, Solv and Cons for a protein
        :param:
        :return
        """

        residue_map = defaultdict(dict)
        max_len = self.length+3
        for r in range(-1,max_len):
            r = str(r)
            residue_map[r] = {'Cum': 0, 'Cum_Dist': 0, 'Cum_Solv': 0, 
			      'Cluster': 0, 'Cluster_Dist': 0, 'Cluster_Solv': 0,
			      'SNAP': 0, 'EVmut': 0, 'Cons': 0, 'Solv': 0}

        # calculate BLOSUM comparison
        self.read_blosum_matrix()
        snap_effect = self.get_snap_effect()
        snap_thresh = self.determine_snap_thresh(snap_effect)
        # print(snap_thresh)
        snap_blosum_diff = self.get_blosum_diff(snap_effect, snap_thresh, 'snap')
        for r in snap_blosum_diff.keys():
            diff = snap_blosum_diff[r]
            r = str(r)
            #print(str(r) + "\t" + str(diff))
            #print(str(r) + "\t" + str(residue_map[r]))
            residue_map[r]['SNAP'] = diff
            #print(str(r) + "\t" + str(residue_map[r]))
        if os.path.isfile(self.evmut_file):
            evc_effect = self.get_evc_effect()
            evc_thresh = self.determine_evc_thresh(evc_effect)
            evc_blosum_diff = self.get_blosum_diff(evc_effect, evc_thresh, 'evmut')
            for r in evc_blosum_diff.keys():
                diff = evc_blosum_diff[r]
                r = str(r)
                residue_map[r]['EVmut'] = diff

        # get per residue conservation
        line_num = 1
        with open(self.cons_file) as cons:
            for line in cons:
                if line_num > 1:
                    splitted = line.split(",")
                    pos = splitted[0]
                    conservation = splitted[2] 
                    residue_map[pos]['Cons'] = conservation
                line_num += 1        

        # get per residue solvent accessibility
        pp = ProfParser(self.solv_file, 'reprof')
        seq = pp.seq
        solv = pp.p_solv_acc
        solv_acc = dict()

        for i in range(0, len(solv)):
            val = solv[i]
            letter = seq[i]
            max_val = self.abs_vals[letter]
            rel_val = val/max_val
            rel_val = round(rel_val, 3)
            pos = i + 1
            residue_map[str(pos)]['Solv'] = rel_val
            solv_acc[pos] = rel_val

        if os.path.isfile(self.evc_file):
            # calculate cumulative coupling scores
            distances = self.get_distances()
            ec_scores = self.get_ec_scores()
            # print(ec_scores)
            self.calc_thresh_avg(ec_scores)
            # print(self.threshold)
            # print(self.average)
            # don't filter
            filtered_scores = self.filter_scores(ec_scores, distances, solv_acc, 'none', 0)
            # print(filtered_scores)
            cumulative_scores = self.calc_cum_scores(filtered_scores)
            # filter by distance
            filtered_dist_scores = self.filter_scores(ec_scores, distances, solv_acc, 'dist', self.cs_dist_thresh)
            cumulative_dist_scores = self.calc_cum_scores(filtered_dist_scores)
            # filter by solvent accessibility
            filtered_solv_scores = self.filter_scores(ec_scores, distances, solv_acc, 'solv', 0)
            cumulative_solv_scores = self.calc_cum_scores(filtered_solv_scores)
        
            # calculate clustering coefficients
            # don't filter
            clustering_coefficients = self.calc_clustering_coefficients(filtered_scores)
            # filter by distance
            filtered_dist_scores = self.filter_scores(ec_scores, distances, solv_acc, 'dist', self.cc_dist_thresh)
            clustering_coefficients_dist = self.calc_clustering_coefficients(filtered_dist_scores)
            # filter by solvent accessibility
            clustering_coefficients_solv = self.calc_clustering_coefficients(filtered_solv_scores)

            for r in range(1,self.length+1):
                cum_score = cumulative_scores[r]
                cum_dist_score = cumulative_dist_scores[r]
                cum_solv_score = cumulative_solv_scores[r]
                cluster_score = clustering_coefficients[r]
                cluster_dist_score = clustering_coefficients_dist[r]
                cluster_solv_score = clustering_coefficients_solv[r]
                r = str(r)
                residue_map[r]['Cum'] = cum_score
                residue_map[r]['Cum_Dist'] = cum_dist_score
                residue_map[r]['Cum_Solv'] = cum_solv_score
                residue_map[r]['Cluster'] = cluster_score
                residue_map[r]['Cluster_Dist'] = cluster_dist_score
                residue_map[r]['Cluster_Solv'] = cluster_solv_score

        #convert to right format
        for r in residue_map.keys():
            r = int(r)
            if r > 0 and r <= self.length:
                r1 = r-2
                r2 = r-1
                r3 = r+1
                r4 = r+2
                r_map = residue_map[str(r)]
                r1_map = residue_map[str(r1)]
                r2_map = residue_map[str(r2)]
                r3_map = residue_map[str(r3)]
                r4_map = residue_map[str(r4)]
                # print(str(r) + "\t" + str(r_map))
                scores = [r_map['Cum'], r2_map['Cum'], r3_map['Cum'], r1_map['Cum'], r4_map['Cum'],
		      r_map['Cum_Solv'], r2_map['Cum_Solv'], r3_map['Cum_Solv'], r1_map['Cum_Solv'], r4_map['Cum_Solv'],
		      r_map['Cum_Dist'], r2_map['Cum_Dist'], r3_map['Cum_Dist'], r1_map['Cum_Dist'], r4_map['Cum_Dist'],
		      r_map['Cluster'], r2_map['Cluster'], r3_map['Cluster'], r1_map['Cluster'], r4_map['Cluster'],
		      r_map['Cluster_Solv'], r2_map['Cluster_Solv'], r3_map['Cluster_Solv'], r1_map['Cluster_Solv'], r4_map['Cluster_Solv'],
		      r_map['Cluster_Dist'], r2_map['Cluster_Dist'], r3_map['Cluster_Dist'], r1_map['Cluster_Dist'], r4_map['Cluster_Dist'],
		      r_map['EVmut'], r2_map['EVmut'], r3_map['EVmut'], r1_map['EVmut'], r4_map['EVmut'],
		      r_map['SNAP'], r2_map['SNAP'], r3_map['SNAP'], r1_map['SNAP'], r4_map['SNAP'],
		      r_map['Solv'], r2_map['Solv'], r3_map['Solv'], r1_map['Solv'], r4_map['Solv'],
		      r_map['Cons'], r2_map['Cons'], r3_map['Cons'], r1_map['Cons'], r4_map['Cons']]
                self.scores[r] = scores


    def get_distances(self):
        line_num = 1
        dist_dict = dict()
        with open(self.dist_file) as f:
            for line in f:
                if line_num > 1:
                    splitted = line.split(",")
                    pos1 = splitted[0]
                    pos2 = splitted[1]
                    comb_pos = pos1 + "/" + pos2
                    dist = 0.0
                    if splitted[len(splitted)-2] != "" and re.search("[a-zA-z]+", splitted[len(splitted)-2]) == 'None':
                        dist = float(splitted[len(splitted)-2])
                    dist_dict[comb_pos] = dist
                line_num += 1
        return dist_dict


    def get_ec_scores(self):
        ec_dict = dict()
        with open(self.evc_file) as f:
            for line in f:
                splitted = line.split(" ")
                comb_pos = splitted[0].strip() + "/" + splitted[2].strip()
                score = float(splitted[4])
                ec_dict[comb_pos] = score
        return ec_dict


    def calc_thresh_avg(self, ec_scores):
        num = self.length_cov*self.factor
        high_scores = list()
        counter = 0
        for key in ec_scores.keys():
            splitted = key.split("/")
            pos1 = int(splitted[0])
            pos2 = int(splitted[1])
            score = ec_scores[key]
            diff = abs(pos1-pos2)
            # print(str(pos1) + "\t" + str(pos2) + "\t" + str(diff))
            if abs(pos1-pos2) > self.seq_dist:
                if len(high_scores) >= num:
                    el = high_scores[0]
                    if el < score:
                        del high_scores[0]
                        high_scores.append(score)
                else:
                    high_scores.append(score)
                high_scores.sort()
            else:
                counter += 1
        #print(counter)
        #print(high_scores)
        #print(len(high_scores))
        thresh = high_scores[0]
        average = round(numpy.mean(high_scores),2)
        self.threshold = thresh
        self.average = average
               

    def filter_scores(self, ec_scores, distances, solv_acc, identifier, dist_thresh):
        filtered_scores = dict()

        for key in ec_scores.keys():
            score = ec_scores[key]
            key_parts = key.split("/")
            pos1 = int(key_parts[0])
            pos2 = int(key_parts[1])
            curr_seq_dist = abs(pos1-pos2)
            if curr_seq_dist > self.seq_dist and score >= self.threshold:
                if identifier == 'none':
                    filtered_scores[key] = score
                elif identifier == 'solv':
                    core = False
                    if pos1 in solv_acc.keys() and pos2 in solv_acc.keys():
                        solv_acc1 = solv_acc[pos1]
                        solv_acc2 = solv_acc[pos2]
                        if solv_acc1 <= 0.1 and solv_acc2 <= 0.1:
                             core = True
                    if core == False:
                        filtered_scores[key] = score
                elif identifier == 'dist':
                    pair_dist = float('Inf')
                    if key in distances:
                        pair_dist = distances[key]
                    if pair_dist <= dist_thresh:
                        filtered_scores[key] = score
        return filtered_scores


    def calc_cum_scores(self, scores):
        """ Calc cumulative coupling scores from a given set of EC scores
        :param scores: list of ec scores
        :return cumulative coupling scores
        """
        cum_scores = dict()
        for i in range(1,self.length+1):
            sum_ec = 0.0
            num = 0
            for j in range(1,self.length+1):
                key = ""
                if i <= j:
                    key = str(i) + "/" + str(j)
                else:
                    key = str(j) + "/" + str(i)
                if key in scores.keys():
                    score = scores[key]
                    sum_ec += score
                    num += 1
            ec_strength = 0.0
            if num > 0:
                ec_strength = sum_ec/self.average
                ec_strength = round(ec_strength,3)
            cum_scores[i] = ec_strength

        return cum_scores        


    def calc_clustering_coefficients(self,scores):
        """ Calc clustering coefficients from a given set of EC scores
        :param scores: list of ec scores
        :return clustering coefficients
        """
        coefficients = dict()
        for i in range(1,self.length+1):
            neighbourhood = set()
            for j in range(1,self.length+1):
                key = ""
                if i <= j:
                    key = str(i) + "/" + str(j)
                else:
                    key = str(j) + "/" + str(i)
                # if key == '228/234':
                #    print(key + "\t" + str(key in scores.keys())) 
                if key in scores.keys():
                    neighbourhood.add(j)

            num_edges = 0
            for n1 in neighbourhood:
                for n2 in neighbourhood:
                    key = ""
                    if n1 <= n2:
                        key = str(n1) + "/" + str(n2)
                    else:
                        key = str(n2) + "/" + str(n1)
                    if key in scores.keys():
                        num_edges += 1
            coeff = 0
            counter = 0
            set_size = len(neighbourhood)
            if set_size > 1:
                coeff = num_edges/(set_size*(set_size-1))
                coeff = round(coeff,3)
                counter =+ 1
            coefficients[i] = coeff
            # print("Clusters > 0: " + str(counter))

        return coefficients


    def get_snap_effect(self):
        """ Read in pre-calculated snap results
        :return: snap scores per residue and mutation
        """
        snap_scores = dict()
        with open(self.snap_file) as f:
            for line in f:
                if "=>" in line:
                    splitted = line.split("=>")
                    identifier = splitted[0].strip()
                    scores = splitted[1].split("\\|")
                    sum_tmp = scores[len(scores)-1].strip()
                    sum_parts = sum_tmp.split("=")
                    sum_final = int(sum_parts[1].strip())
                    snap_scores[identifier] = sum_final
                    # print(identifier + "\t" + str(sum_final))
        return snap_scores


    def determine_snap_thresh(self, scores):
        """Determine the smallest value for SNAP2 still to be classified as having an effect
        :param scores
        :return smallest value
        """
        values = list()
        for key in scores.keys():
            values.append(scores[key])
        sorted_values = sorted(values)
        index = int(len(values)*(1-self.snap_percentage))
        thresh = sorted_values[index]

        return thresh


    def get_evc_effect(self):
        """ Read in pre-calculated EVmutation results
        :return EVmutation scores per residue and mutation
        """
        evc_scores = dict()
        with open(self.evmut_file) as f:
            for line in f:
                splitted = line.strip().split()
                #print(splitted)
                key = splitted[0]
                #print(key)
                value = float(splitted[1])
                evc_scores[key] = value
        return evc_scores


    def determine_evc_thresh(self, scores):
        """Determine the smallest value for EVmutation still to be classified as having an effect
        :param scores
        :return smallest value
        """
        positive_scores = list()
        for key in scores.keys():
            val = scores[key]
            val = abs(val)
            positive_scores.append(val)
        sorted_scores = sorted(positive_scores)
        size = len(sorted_scores)
        num_scores = int(size*self.evc_percentage)
        index = size - num_scores
        thresh = sorted_scores[index]
        return thresh


    def get_blosum_diff(self, scores, thresh, method):
        """ Calculate difference between BLOSUM and effect predictions
        :param scores: effect prediction scores
        :param thresh: smallest value to still classify as having an effect
        :param method: method used to predict effects
        :return BLOSUM difference for each residue
        """
        blosum_diff = dict()
        for i in range(1,self.length+1):
            prog = re.compile("[A-Z]"+str(i)+"[A-Z]")
            num_accepted = 0
            num_mutations = 0
            for mut in scores.keys():
                #print(mut)
                res = prog.search(mut)
                #print(res)
                if res != None:
                    aa1 = mut[:1]
                    aa2 = mut[len(mut)-1:]
                    blosum_score = self.get_blosum_score(aa1, aa2)
                    #print(mut + "\t" + aa1 + "\t" + aa2 + "\t" + str(blosum_score))

                    if blosum_score >= self.blosum_thresh:
                        num_accepted += 1
                        mut_score = scores[mut]
                        if method == 'evmut':
                            mut_score = abs(mut_score)
                        if mut_score >= thresh:
                            num_mutations += 1
            diff = 0.0
            if num_accepted > 0:
                diff = num_mutations/num_accepted
                diff = round(diff, 3)
            blosum_diff[i] = diff
            #print(str(i) + "\t" + str(diff))
        return blosum_diff


    def read_blosum_matrix(self):
        line_num = 1
        pos_mapping = dict()
        with open(self.blosum_file) as f:
            for line in f:
                splitted = line.strip().split()
                if line_num == 1:
                    for i in range(0,len(splitted)):
                        pos = i + 1
                        pos_mapping[pos] = splitted[i]
                else:
                    pos1 = splitted[0]
                    for j in range(1,len(splitted)):
                        pos2 = pos_mapping[j]
                        self.blosum_mat[pos1][pos2] = int(splitted[j])
                        # print(pos1 + "\t" + pos2 + "\t" + splitted[j])
                line_num += 1


    def get_blosum_score(self, aa1, aa2):
        return self.blosum_mat[aa1][aa2]


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
        
