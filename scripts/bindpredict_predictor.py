import os
import glob
import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn import model_selection
from sklearn.model_selection import PredefinedSplit
from sklearn import metrics
from sklearn.externals import joblib
from imblearn.over_sampling import SMOTE
from scripts.protein import Protein
from scripts.file_manager import FileManager
# from scripts.reader import Reader
from scripts.select_model import SelectModel
from scripts.helper import Helper

class Predictor(object):    
    def __init__(self):
        self.query_proteins = dict()
        self.helper = Helper()
        self.f = FileManager()

    def run_prediction(self,args):
        # 1) Read files
        self.get_snap_files(args.snap_folder,args.snap_suffix)
        self.get_evc_files(args.evc_folder)
        self.get_solv_files(args.solv_acc_folder,args.profacc_suffix)
        # 2) Load saved model
        cols_to_remove, ranges, words = SelectModel.define_model(args.model)
        classifier = joblib.load("model/" + args.model + ".pkl")
        # 3) Prepare predictions
        self.calculate_scores(self.query_proteins)    
        for p in self.query_proteins:
            protein_scores = self.prepare_scores_prediction(p, cols_to_remove, ranges, words)
            # 4) Run predictions
            predictions = classifier.predict_proba(protein_scores)
            # 5) Write output
            self.write_output_predictions(args.output_folder, p, predictions)


    def train_predictor(self,args):
        run = int(args.run)
        # 1) Read files
        print('Read files...')
        self.get_snap_files(args.snap_folder,args.snap_suffix)
        self.get_evc_files(args.evc_folder)
        self.get_solv_files(args.solv_acc_folder,args.profacc_suffix)
        self.get_binding_site_information()
        # 2) Get all needed scores per protein and per residue
        print('Calculate scores...')
        ids_folder = self.f.ids_folder
        ids_file_prefix = self.f.ids_file
        split1_ids, split2_ids, split3_ids, split4_ids, split5_ids = self.helper.read_id_lists(ids_folder, ids_file_prefix) 
        # print(len(split1_ids))
        self.calculate_scores(split1_ids)
        self.calculate_scores(split2_ids)
        self.calculate_scores(split3_ids)
        self.calculate_scores(split4_ids)
        self.calculate_scores(split5_ids)
        # 3) Determine model to develop
        print('Prepare model, data and splits...')
        cols_to_remove, ranges, words = SelectModel.define_model(args.model)
        # 4) Prepare predictions
        x_split1, y_split1, length_split1 = self.prepare_splits(split1_ids, cols_to_remove, ranges, words)
        x_split2, y_split2, length_split2 = self.prepare_splits(split2_ids, cols_to_remove, ranges, words)
        x_split3, y_split3, length_split3 = self.prepare_splits(split3_ids, cols_to_remove, ranges, words)
        x_split4, y_split4, length_split4 = self.prepare_splits(split4_ids, cols_to_remove, ranges, words)
        x_split5, y_split5, length_split5 = self.prepare_splits(split5_ids, cols_to_remove, ranges, words)
        # 5) Determine train-cross-train-test-split
        if run == 1:
            x_train = x_split1 + x_split2 + x_split3
            y_train = y_split1 + y_split2 + y_split3
            x_cross = x_split4
            y_cross = y_split4
            x_test = x_split5
            y_test = y_split5
            length_test = length_split5
        if run == 2:
            x_train = x_split5 + x_split1 + x_split2
            y_train = y_split5 + y_split1 + y_split2
            x_cross = x_split3
            y_cross = y_split3
            x_test = x_split4
            y_test = y_split4
            length_test = length_split4
        if run == 3:
            x_train = x_split4 + x_split5 + x_split1
            y_train = y_split4 + y_split5 + y_split1
            x_cross = x_split2
            y_cross = y_split2
            x_test = x_split3
            y_test = y_split3
            length_test = length_split3
        if run == 4:
            x_train = x_split3 + x_split4 + x_split5
            y_train = y_split3 + y_split4 + y_split5
            x_cross = x_split1
            y_cross = y_split1
            x_test = x_split2
            y_test = y_split2
            length_test = length_split2
        if run == 5:
            x_train = x_split2 + x_split3 + x_split4
            y_train = y_split2 + y_split3 + y_split4
            x_cross = x_split5
            y_cross = y_split5
            x_test = x_split1
            y_test = y_split1
            length_test = length_split1
        # 6) Load more data if specified
        if words[1] == "yes":
            print('Load more data...')
            ids_more = self.helper.read_more_data(ids_folder, run)
            ids_more = ids_more[0:578]	#ensure that amount of more data is the same size for all runs
            self.calculate_scores(ids_more)
            x_more, y_more, length_more = self.prepare_splits(ids_more, cols_to_remove, ranges, words)
            x_train = x_train + x_more
            y_train = y_train + y_more
        # 7) Train model
        print('Train model...')
        for el in x_train:
            del el[len(el)-1]

        # print(x_cross)
        for el in x_cross:
            del el[len(el)-1]

        # print(x_train)
        sm = SMOTE(random_state=42)
        x_train, y_train = sm.fit_sample(x_train, y_train)
        x_train = np.array(x_train)
        x_cross = np.array(x_cross)
        y_train = np.array(y_train)
        y_cross = np.array(y_cross)

        x_whole = np.concatenate((x_train, x_cross), axis=0)
        y_whole = np.concatenate((y_train, y_cross), axis=0)

        test_fold_array = [-1]*len(y_train) + [0]*len(y_cross)
        ps = PredefinedSplit(test_fold_array)
        params = {'alpha': [0.0001], 'random_state': [1], 'hidden_layer_sizes': [(100,)], 'solver': ['sgd'], 'learning_rate_init': [0.00002],
		  'momentum': [0.018], 'max_iter': [500]}
        grid = model_selection.GridSearchCV(MLPClassifier(), params, cv=ps)
        grid.fit(x_whole, y_whole)
        classifier = grid.best_estimator_
        # 8) Calculate predictions and get performance
        print("Calculate predictions and performance...")
        test_res = []
        for index in range(0,len(x_test)):
            el = x_test[index]
            print(el)
            label = y_test[index]
            to_test = el[:len(el)-1]
            # set distance values = 0, TODO immediately remove afterwards!!
            pos = el[len(el)-1]
            to_test = np.reshape(to_test, (1,-1))
            proba = classifier.predict_proba(to_test)
            prediction = proba[0][1]
            classify = ""
            if prediction < 0.6:
                if label == 0:
                    classify = "TN"
                else:
                    classify = "FN"
            else:
                if label == 0:
                    classify = "FP"
                else:
                    classify = "TP"
            res = [pos, prediction, classify]
            test_res.append(res) 
        self.calculate_performance(test_res)
        # 9) Write output
        print('Write output...')
        self.write_output_cv(args.output, run, test_res)
        self.write_output_proba(args.output, run, test_res)


    def save_trained_model(self,args):
        # 1) Read files
        self.get_snap_files(args.snap_folder,args.snap_suffix)
        self.get_evc_files(args.evc_folder)
        self.get_solv_files(args.solv_folder,args.profacc_suffix)
        self.get_binding_site_information()
        # 2) Determine model to train and save
        cols_to_remove, ranges, words = SelectModel.define_model(args.model)
        # 3) Train model using all splits (+ more data?)
        ids_folder = self.f.ids_folder
        ids_file_prefix = self.f.ids_file
        split1_ids, split2_ids, split3_ids, split4_ids, split5_ids = self.helper.read_id_lists(ids_folder, ids_file_prefix) 
        self.calculate_scores(split1_ids)
        self.calculate_scores(split2_ids)
        self.calculate_scores(split3_ids)
        self.calculate_scores(split4_ids)
        self.calculate_scores(split5_ids)
        
        x_split1, y_split1, length_split1 = self.prepare_splits(split1_ids, cols_to_remove, ranges, words)
        x_split2, y_split2, length_split2 = self.prepare_splits(split2_ids, cols_to_remove, ranges, words)
        x_split3, y_split3, length_split3 = self.prepare_splits(split3_ids, cols_to_remove, ranges, words)
        x_split4, y_split4, length_split4 = self.prepare_splits(split4_ids, cols_to_remove, ranges, words)
        x_split5, y_split5, length_split5 = self.prepare_splits(split5_ids, cols_to_remove, ranges, words)
        x_data = x_split1 + x_split2 + x_split3 + x_split4 + x_split5
        y_data = y_split1 + y_split2 + y_split3 + y_split4 + y_split5
        # 4) Load more data if specified
        if words[1] == "yes":
            ids_more1 = self.helper.read_more_data(ids_folder, 1)
            ids_more1 = ids_more1[0:578]	#ensure that amount of more data is the same size for all runs
            self.calculate_scores(ids_more1)
            x_more1, y_more1, length_more1 = self.prepare_splits(ids_more1, cols_to_remove, ranges, words)
            ids_more2 = self.helper.read_more_data(ids_folder, 1)
            ids_more2 = ids_more1[0:578]	#ensure that amount of more data is the same size for all runs
            self.calculate_scores(ids_more2)
            x_more2, y_more2, length_more2 = self.prepare_splits(ids_more2, cols_to_remove, ranges, words)
            ids_more3 = self.helper.read_more_data(ids_folder, 1)
            ids_more3 = ids_more1[0:578]	#ensure that amount of more data is the same size for all runs
            self.calculate_scores(ids_more3)
            x_more3, y_more3, length_more3 = self.prepare_splits(ids_more3, cols_to_remove, ranges, words)
            ids_more4 = self.helper.read_more_data(ids_folder, 1)
            ids_more4 = ids_more1[0:578]	#ensure that amount of more data is the same size for all runs
            self.calculate_scores(ids_more4)
            x_more4, y_more4, length_more4 = self.prepare_splits(ids_more4, cols_to_remove, ranges, words)
            ids_more5 = self.helper.read_more_data(ids_folder, 1)
            ids_more5 = ids_more1[0:578]	#ensure that amount of more data is the same size for all runs
            self.calculate_scores(ids_more5)
            x_more5, y_more5, length_more5 = self.prepare_splits(ids_more5, cols_to_remove, ranges, words)
            x_data = x_data + x_more1 + x_more2 + x_more3 + x_more4 + x_more5
            y_data = y_data + y_more1 + y_more2 + y_more3 + y_more4 + y_more5
        # 5) Train model
        for el in x_data:
            del el[len(el)-1]

        sm = SMOTE(random_state=42)
        x_data, y_data = sm.fit_sample(x_data, y_data)
        classifier = MLPClassifier(alpha=0.0001, random_state=1, hidden_layer_sizes=[(100,)], solver='sgd', learning_rate_init=0.00002, momentum=0.018, max_iter=500, early_stopping=True)
        classifier.fit(x_data, y_data)
        # 6) Save trained model 
        model_path = "model/" + args.model + ".pkl" 
        joblib.dump(classifier, model_path)


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
            name = name[:suffix_start]	# reduce the file name to the prefix only
            self.query_proteins[name] = Protein()
            self.query_proteins[name].snap_file = f
            #print(self.f.blosum_file())
            self.query_proteins[name].blosum_file = self.f.blosum_file()            


    def get_evc_files(self, folder):
        """ For each id, determine the different evc files.
        :param folder: Folder where the folders with EVcouplings predictions are stored
        :return:
        """

        for p in self.query_proteins:
            evc_folder = os.path.join(folder, p)
            evc_file = os.path.join(evc_folder, (p + self.f.ec_suffix))
            dist_file = os.path.join(evc_folder, (p + self.f.dist_suffix))
            evmut_file = os.path.join(evc_folder, (p + self.f.evmut_suffix))
            cons_file = os.path.join(evc_folder, (p + self.f.cons_suffix))
            align_file = os.path.join(evc_folder, (p + self.f.align_suffix))
            
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
        :param folder: Folder where PROFacc predictions are stored
        :param suffix: suffix of the files (default: reprof)
        :return:
        """

        solv_pattern = os.path.join(folder, suffix)
        for f in glob.iglob(solv_pattern):
            name = os.path.basename(f)
            suffix_start = name.find('.')
            name = name[:suffix_start]
            self.query_proteins[name].solv_file = f


    def get_binding_site_information(self):
        """ Read in the binding site information for all proteins.
        :param
        :return
        """

        bind_res_file = self.f.binding_file
        with open(bind_res_file) as f:
            for row in f:
                splitted = row.strip().split("\t")
                uni_id = splitted[0]
                residues = splitted[1:]
                if uni_id in self.query_proteins:
                    self.query_proteins[uni_id].binding_res = residues
                    # print(uni_id)
                    # print(self.query_proteins[uni_id].binding_res)
                

    def calculate_scores(self, id_list):
        """ Calculate cum scores, clustering coefficients, SNAP2/BLOSUM, EVmut/BLOSUM, Solv and Cons.
        :param id_list: list of ids to calculate scores for
        :return
        """
        
        for i in id_list:
            self.query_proteins[i].calc_scores()


    def prepare_splits(self, id_list, cols_to_remove, ranges, words):
        """For each id, get scores needed for the given model and get corrrect labelling for each residue
        :param id_list: list of ids
        :param cols_to_remove: scores not to use for the model
        :param ranges: scores to calculate average for if model is an average model
        :param words: specifying whether more data should be included and an average model is used 
        :return data: scores for each residue
        :return labels: true label for each residue
        :return lengths: length of each protein
        """

        data = []
        labels = []
        lengths = []

        for i in id_list:
            scores = self.query_proteins[i].scores
            binding = self.query_proteins[i].binding_res
            for res in scores.keys():
                resi = str(res)
                scores_res = scores[res]            
                scores_res_new = [float(a) for j,a in enumerate(scores_res) if j not in cols_to_remove]
                i_res = i + "_" + str(res)
                scores_res_new = scores_res_new + [i_res] 
                data.append(scores_res_new)
                if resi in binding:	# residue is annotated as binding
                    labels.append(1)
                else:
                    labels.append(0)
            length = self.query_proteins[i].length
            lengths.append(length)

        return data, labels, lengths


    def prepare_scores_prediction(self, p, cols_to_remove, ranges, words):
        data = []
        scores = self.query_proteins[p].scores
        for res in scores.keys():
            scores_res = scores[res]
            row = [i for j, i in enumerate(scores_res) if j not in cols_to_remove]
            data.append(row)

        return data
    

    def calculate_performance(self, test_results):
        """ Calculate common performance measurements for all proteins
        :param test_results: Prediction results for all residues
        :return
        """

        id_set = set()
        for t in test_results:
            splitted = t[0].split("_")
            protein = splitted[0]
            residue = splitted[1]
            label = t[2]
            self.query_proteins[protein].performance[label] += 1
            id_set.add(protein)

        id_list = list(id_set)
        for i in id_list:
            tp = self.query_proteins[i].performance['TP']
            fp = self.query_proteins[i].performance['FP']
            tn = self.query_proteins[i].performance['TN']
            fn = self.query_proteins[i].performance['FN'] 
            cov = 0
            if tp > 0 or fn > 0:
                cov = tp/(tp+fn)
                cov = round(cov, 3)
            prec = 0
            if tp > 0 or fp > 0:
                prec = tp/(tp+fp)
                prec = round(prec, 3)
            acc = (tp+tn)/(tp+fp+tn+fn)
            acc = round(acc, 3)
            f1 = 0
            if cov > 0 or prec > 0:
                f1 = 2*cov*prec/(cov+prec)
                f1 = round(f1, 3)

            self.query_proteins[i].performance['Cov'] = cov
            self.query_proteins[i].performance['Prec'] = prec
            self.query_proteins[i].performance['Acc'] = acc
            self.query_proteins[i].performance['F1'] = f1

    
    def write_output_predictions(self, folder, p, predictions):
        """ Write predictions results into output file.
        :param folder: Folder to write the predictions results to
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
            
    def write_output_proba(self, out_file, run, results):
        """ Write probabilities results into output file.
        :param file: File to write results to
        :param run: Run number to use as suffix
        :return
        """

        with open((out_file + "." + str(run) + ".probabilities"), 'w') as out:
            for r in results:
                # format: pos\tproba\tlabel
                out.write(r[0] + "\t" + str(r[1]) + "\t" + r[2] + "\n")


    def write_output_cv(self, out_file, run, results):
        """ Write performance results into output file.
        :param file: File to write the performance results to
        :return:
        """
        
        proteins = set()
        for r in results:
            p = r[0].split('_')[0]
            proteins.add(p)

        proteins = list(proteins)

        with open((out_file + "." + str(run) + ".performance"), 'w') as out:
            out.write("ID\tTP\tFP\tTN\tFN\tCov\tPrec\tF1\tAcc\n")
            for p in proteins:
                performance = self.query_proteins[p].performance
                out.write(p + "\t" + str(performance['TP']) + "\t" + str(performance['FP']) + "\t" + str(performance['TN']) + "\t" + str(performance['FN']) +
                          "\t" + str(performance['Cov']) + "\t" + str(performance['Prec']) + "\t" + str(performance['F1']) + "\t" + str(performance['Acc']) + "\n")
                 
 
def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)
