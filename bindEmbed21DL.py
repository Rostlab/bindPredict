from data_preparation import ProteinInformation
from ml_trainer import MLTrainer
from ml_predictor import MLPredictor
from config import FileSetter, FileManager

import numpy as np
from sklearn.model_selection import PredefinedSplit


class BindEmbed21DL(object):

    @staticmethod
    def cross_train_pipeline(params, model_output, predictions_output, ri):
        """
        Run cross-training pipeline for a specific set of parameters
        :param params:
        :param model_output: If None, trained model is not written
        :param predictions_output: If None, predictions are not written
        :param ri: Should RI or raw probabilities be written?
        :return:
        """

        print("Prepare data")
        ids = []
        fold_array = []
        for s in range(1, 6):
            ids_in = '{}{}.txt'.format(FileSetter.split_ids_in(), s)
            split_ids = FileManager.read_ids(ids_in)

            ids += split_ids
            fold_array += [s] * len(split_ids)

        ps = PredefinedSplit(fold_array)

        # get sequences + maximum length + labels
        sequences, max_length, labels = ProteinInformation.get_data(ids)
        embeddings = FileManager.read_embeddings(FileSetter.embeddings_input())

        proteins = dict()
        trainer = MLTrainer(pos_weights=params['weights'])

        for train_index, test_index in ps.split():

            split_counter = fold_array[test_index[0]]
            train_ids = [ids[train_idx] for train_idx in train_index]
            validation_ids = [ids[test_idx] for test_idx in test_index]

            print("Train model")
            model_split = trainer.train_validate(params, train_ids, validation_ids, sequences, embeddings, labels,
                                                 max_length, verbose=False)

            if model_output is not None:
                model_path = '{}{}.pt'.format(model_output, split_counter)
                FileManager.save_classifier_torch(model_split, model_path)

            print("Calculate predictions per protein")
            ml_predictor = MLPredictor(model_split)
            curr_proteins = ml_predictor.predict_per_protein(validation_ids, sequences, embeddings, labels, max_length)

            proteins = {**proteins, **curr_proteins}

            if predictions_output is not None:
                FileManager.write_predictions(proteins, predictions_output, 0.5, ri)

        return proteins

    @staticmethod
    def hyperparameter_optimization_pipeline(params, num_splits, result_file):
        """
        Development pipeline used to optimize hyperparameters
        :param params:
        :param num_splits:
        :param result_file:
        :return:
        """

        print("Prepare data")
        ids = []
        fold_array = []
        for s in range(1, num_splits + 1):
            ids_in = '{}{}.txt'.format(FileSetter.split_ids_in(), s)
            split_ids = FileManager.read_ids(ids_in)

            ids += split_ids
            fold_array += [s] * len(split_ids)

        ids = np.array(ids)

        # get sequences + maximum length + labels
        sequences, max_length, labels = ProteinInformation.get_data(ids)
        embeddings = FileManager.read_embeddings(FileSetter.embeddings_input())

        print("Perform hyperparameter optimization")
        trainer = MLTrainer(pos_weights=params['weights'])
        del params['weights']  # remove weights to not consider as parameter for optimization

        model = trainer.cross_validate(params, ids, fold_array, sequences, embeddings, labels, max_length, result_file)

        return model

    @staticmethod
    def prediction_pipeline(model_prefix, cutoff, result_folder, ids, sequences, ri):
        """
        Run predictions with bindEmbed21DL for a given list of proteins
        :param model_prefix:
        :param cutoff: Cutoff to use to define prediction as binding (default: 0.5)
        :param result_folder:
        :param ids:
        :param sequences:
        :param ri: Should RI or raw probabilities be written?
        :return:
        """

        print("Prepare data")
        embeddings = FileManager.read_embeddings(FileSetter.embeddings_input())

        proteins = dict()
        for i in range(0, 5):
            print("Load model")
            model_path = '{}{}.pt'.format(model_prefix, i + 1)
            model = FileManager.load_classifier_torch(model_path)

            print("Calculate predictions")
            ml_predictor = MLPredictor(model)
            curr_proteins = ml_predictor.predict_per_protein(ids, sequences, embeddings, labels, max_length)

            for k in curr_proteins.keys():
                if k in proteins.keys():
                    prot = proteins[k]
                    prot.add_predictions(curr_proteins[k].predictions)
                else:
                    proteins[k] = curr_proteins[k]

        for k in proteins.keys():
            proteins[k].normalize_predictions(5)

        if result_folder is not None:
            FileManager.write_predictions(proteins, result_folder, cutoff, ri)

        return proteins
