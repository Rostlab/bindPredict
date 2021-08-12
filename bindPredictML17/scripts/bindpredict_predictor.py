import os
import sys
import glob
import numpy as np
from sklearn.externals import joblib
from scripts.protein import Protein
from scripts.file_manager import FileManager
from scripts.select_model import SelectModel


class Predictor(object):
    def __init__(self):
        self.query_proteins = dict()
        self.fm = FileManager()

    def run_prediction(self, args):
        # 1) Read files
        print('1) Read files...')
        self.get_snap_files(args.snap_folder, args.snap_suffix)
        self.get_evc_files(args.evc_folder)
        self.get_solv_files(args.solv_acc_folder, args.profacc_suffix)
        # 2) Load saved model
        print('2) Load saved model...')
        cols_to_remove, ranges, words = SelectModel.define_model(args.model)
        classifier = joblib.load("trained_model/" + args.model + ".pkl")
        # 3) Prepare predictions
        print('3) Calculate scores for query...')
        self.calculate_scores(self.query_proteins)
        print('4) Run predictions...')
        for p in self.query_proteins:
            protein_scores = self.prepare_scores_prediction(p, cols_to_remove)
            protein_scores = np.array(protein_scores)
            # 4) Run predictions per protein
            predictions = []
            for index in range(0, len(protein_scores)):
                el = protein_scores[index]
                to_test = el[:len(el) - 1]
                pos = el[len(el) - 1]
                to_test = np.reshape(to_test, (1, -1))
                to_test = to_test.astype(float)
                proba = classifier.predict_proba(to_test)
                prediction = proba[0][1]
                row = [pos, prediction]
                predictions = predictions + [row]
            # 5) Write output
            self.write_output_predictions(args.output_folder, predictions, p)

    def get_snap_files(self, folder, suffix):
        """ Read all ids in a big dictionray.
        :param folder: Folder where the snap2 predictions are stored
        :param suffix: suffix of the files (default: snap2)
        :return:
        """

        snap_pattern = os.path.join(folder, suffix)
        for f in glob.iglob(snap_pattern):
            name = os.path.basename(f)
            suffix_start = name.find('.')
            name = name[:suffix_start]  # reduce the file name to the prefix only
            self.query_proteins[name] = Protein()
            self.query_proteins[name].snap_file = f
            self.query_proteins[name].blosum_file = self.fm.blosum_file

    def get_evc_files(self, folder):
        """ For each id, determine the different evc files.
        :param folder: Folder where the folders with EVcouplings predictions are stored
        :return:
        """

        for p in self.query_proteins:
            evc_folder = os.path.join(folder, p)
            evc_file = os.path.join(evc_folder, (p + self.fm.ec_suffix))
            dist_file = os.path.join(evc_folder, (p + self.fm.dist_suffix))
            evmut_file = os.path.join(evc_folder, (p + self.fm.evmut_suffix))
            cons_file = os.path.join(evc_folder, (p + self.fm.cons_suffix))
            align_file = os.path.join(evc_folder, (p + self.fm.align_suffix))

            self.query_proteins[p].evc_file = evc_file
            self.query_proteins[p].dist_file = dist_file
            self.query_proteins[p].evmut_file = evmut_file
            self.query_proteins[p].cons_file = cons_file

            counter = 1
            with open(align_file) as af:
                for line in af:
                    if counter == 2:
                        splitted = line.split(",")
                        length = int(splitted[3])
                        length_cov = int(splitted[4])
                        self.query_proteins[p].length = length
                        self.query_proteins[p].length_cov = length_cov
                    counter += 1

    def get_solv_files(self, folder, suffix):
        """ For each id, safe the solv file.
        :param folder: Folder where PROFacc/Reprof/PROFphd predictions are stored
        :param suffix: suffix of the files (default: reprof)
        :return:
        """

        solv_pattern = os.path.join(folder, suffix)
        for f in glob.iglob(solv_pattern):
            name = os.path.basename(f)
            suffix_start = name.find('.')
            name = name[:suffix_start]
            self.query_proteins[name].solv_file = f

    def calculate_scores(self, id_list):
        """ Calculate cum scores, clustering coefficients, SNAP2/BLOSUM, EVmut/BLOSUM, Solv and Cons.
        :param id_list: list of ids to calculate scores for
        :return
        """

        for i in id_list:
            self.query_proteins[i].calc_scores()

    def prepare_scores_prediction(self, p, cols_to_remove):
        data = []
        scores = self.query_proteins[p].scores
        for res in scores.keys():
            resi = str(res)
            scores_res = scores[res]
            row = [float(i) for j, i in enumerate(scores_res) if j not in cols_to_remove]
            row = row + [resi]
            data.append(row)

        return data

    @staticmethod
    def write_output_predictions(folder, predictions, p):
        """ Write predictions results into output file.
        :param folder: Folder to write the prediction results to
        :param p: Protein id to use as filename
        :param predictions: Per residue predictions
        :return:
        """

        out_file = os.path.join(folder, (p + ".bindPredict_out"))
        with open(out_file, 'w') as out:
            # format: res, proba_pos
            for pred in predictions:
                label = "nb"
                if pred[1] >= 0.6:
                    label = "b"
                out.write(str(pred[0]) + "\t" + str(pred[1]) + "\t" + label + "\n")


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
