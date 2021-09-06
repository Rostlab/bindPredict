from bindEmbed21DL import BindEmbed21DL
from assess_performance import PerformanceAssessment
from config import FileSetter, FileManager, GeneralInformation
from data_preparation import ProteinInformation

import sys
from pathlib import Path


def main():
    GeneralInformation.seed_all(42)

    keyword = sys.argv[1]

    path = ''  # TODO set path to working directory
    Path(path).mkdir(parents=True, exist_ok=True)

    cutoff = 0.5
    ri = False  # Should RI or raw probabilities be written?

    sequences = FileManager.read_fasta(FileSetter.fasta_file())

    if keyword == 'optimize-architecture':
        params = {'lr': [0.1, 0.05, 0.01], 'betas': [(0.9, 0.999)], 'eps': [1e-8], 'weight_decay': [0.0],
                  'features': [512, 256, 128, 64, 32], 'kernel': [3, 5, 7], 'stride': [1],
                  'dropout': [0.7, 0.5, 0.4, 0.3, 0.2], 'epochs': [200], 'early_stopping': [True],
                  'weights': [[8.9, 7.7, 4.4]]}

        result_file = '{}/cross_validation_results.txt'.format(path)

        BindEmbed21DL.hyperparameter_optimization_pipeline(params, 5, result_file)

    elif keyword == 'best-training':
        params = {'lr': 0.01, 'betas': (0.9, 0.999), 'eps': 1e-8, 'weight_decay': 0.0, 'features': 128, 'kernel': 5,
                  'stride': 1, 'dropout': 0.7, 'epochs': 200, 'early_stopping': True, 'weights': [8.9, 7.7, 4.4]}

        model_prefix = '{}/trained_model'.format(path)
        prediction_folder = '{}/predictions'.format(path)
        Path(prediction_folder).mkdir(parents=True, exist_ok=True)

        proteins = BindEmbed21DL.cross_train_pipeline(params, model_prefix, prediction_folder, ri)

        # assess performance
        # labels = ProteinInformation.get_labels(proteins.keys(), sequences)
        # model_performances = PerformanceAssessment.combine_protein_performance(proteins, cutoff, labels)
        # PerformanceAssessment.print_performance_results(model_performances)

    elif keyword == 'testing':
        model_prefix = '{}/trained_model'.format(path)
        prediction_folder = '{}/predictions_testset/'.format(path)
        Path(prediction_folder).mkdir(parents=True, exist_ok=True)

        ids_in = FileSetter.test_ids_in()
        fasta_file = FileSetter.fasta_file()
        proteins = BindEmbed21DL.prediction_pipeline(model_prefix, cutoff, prediction_folder, ids_in, fasta_file, ri)

        # assess performance
        # labels = ProteinInformation.get_labels(proteins.keys(), sequences)
        # model_performances = PerformanceAssessment.combine_protein_performance(proteins, cutoff, labels)
        # PerformanceAssessment.print_performance_results(model_performances)


main()
