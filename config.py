from pandas import DataFrame
import torch
import os
import numpy
import random
from collections import defaultdict


class FileSetter(object):

    @staticmethod
    def t5_dir():
        return ''  # TODO set path to embeddings, embeddings should be in .npy-format with one embedding per

    @staticmethod
    def predictions_folder():
        return ''  # TODO set path to where predictions should be written

    @staticmethod
    def profile_db():
        # TODO set path to pre-computed big_80 database
        # can be downloaded from ftp://rostlab.org/bindEmbed21/profile_db.tar.gz
        return ''

    @staticmethod
    def lookup_fasta():
        # TODO set path to FASTA file of lookup set
        # can be downloaded from ftp://rostlab.org/bindEmbed21/lookup.fasta
        return ''

    @staticmethod
    def lookup_db():
        # TODO set path to pre-computed lookup database
        # can be downloaded from ftp://rostlab.org/bindEmbed21/lookup_db.tar.gz
        return ''

    @staticmethod
    def mmseqs_output():
        # TODO set path to where MMseqs2 output folder should be written.
        # tmp files will be stored in this folder in a sub-directory tmp/
        # predictions will be stored in this folder in a sub-directory hbi_predictions/
        return ''

    @staticmethod
    def query_set():
        return ''  # TODO set path to FASTA set of query sequences to generate predictions for

    @staticmethod
    def mmseqs_path():
        return ''  # TODO set path to MMseqs2 installation

    @staticmethod
    def split_ids_in():
        return 'data/development_set/ids_split'  # cv splits used during development; available on GitHub

    @staticmethod
    def test_ids_in():
        return 'data/development_set/uniprot_test.txt'  # test ids used during development; available on GitHub

    @staticmethod
    def fasta_file():
        return 'data/development_set/all.fasta'  # path to development set; available on GitHub

    @staticmethod
    def binding_residues_by_ligand(ligand):
        return 'data/development_set/binding_residues_2.5_{}.txt'.format(ligand)
        # files with binding labels used during development; available on GitHub


class FileManager(object):

    @staticmethod
    def read_ids(file_in):
        """
        Read list of ids into list
        :param file_in:
        :return:
        """
        ids = []
        with open(file_in) as read_in:
            for line in read_in:
                ids.append(line.strip())

        return ids

    @staticmethod
    def read_fasta(file_in):
        """
        Read sequences from FASTA file
        :param file_in:
        :return: dict with key: ID, value: sequence
        """
        sequences = dict()
        current_id = None

        with open(file_in) as read_in:
            for line in read_in:
                line = line.strip()
                if line.startswith(">"):
                    current_id = line[1:]
                    sequences[current_id] = ''
                else:
                    sequences[current_id] += line

        return sequences

    @staticmethod
    def read_binding_residues(file_in):
        """
        Read binding residues from file
        :param file_in:
        :return:
        """
        binding = dict()

        with open(file_in) as read_in:
            for line in read_in:
                splitted_line = line.strip().split()
                if len(splitted_line) > 1:
                    identifier = splitted_line[0]
                    residues = splitted_line[1].split(',')
                    residues_int = [int(r) for r in residues]

                    binding[identifier] = residues_int

        return binding

    @staticmethod
    def read_mmseqs_alignments(file_in):
        """Read MMseqs2 alignments"""

        mmseqs = defaultdict(defaultdict)

        with open(file_in) as read_in:
            for line in read_in:
                splitted_line = line.strip().split()
                query_id = splitted_line[0]
                target_id = splitted_line[1]
                qstart = int(splitted_line[5])
                tstart = int(splitted_line[6])
                qaln = splitted_line[7]
                taln = splitted_line[8]

                mmseqs[query_id][target_id] = {'qstart': qstart, 'tstart': tstart, 'qaln': qaln, 'taln': taln}

        return mmseqs

    @staticmethod
    def save_cv_results(cv_results, file):
        """
        Save CV results to csv file
        :param cv_results:
        :param file:
        :return:
        """

        cv_dataframe = DataFrame.from_dict(cv_results)
        cv_dataframe.to_csv(path_or_buf=file)

    @staticmethod
    def save_classifier_torch(classifier, model_path):
        """Save pre-trained model"""
        torch.save(classifier, model_path)

    @staticmethod
    def load_classifier_torch(model_path):
        """ Load pre-saved model """
        if torch.cuda.is_available():
            device = 'cuda:0'
        else:
            device = 'cpu'
        classifier = torch.load(model_path, map_location=device)
        return classifier

    @staticmethod
    def write_predictions(proteins, out_folder, cutoff, ri):
        """
        Write predictions for a set of proteins
        :param proteins:
        :param out_folder:
        :param cutoff: Cutoff to define whether a residue is binding or not
        :param ri: Should raw probabilities or RI be written to file?
        :return:
        """
        for k in proteins.keys():
            p = proteins[k]
            predictions = p.predictions
            predictions = predictions.squeeze()
            out_file = os.path.join(out_folder, (k + '.bindPredict_out'))

            FileManager.write_predictions_single_protein(out_file, predictions, cutoff, ri)

    @staticmethod
    def write_predictions_single_protein(out_file, predictions, cutoff, ri):
        """ Write predictions for a specific protein """
        with open(out_file, 'w') as out:
            if ri:
                out.write("Position\tMetal.RI\tMetal.Class\tNuc.RI\tNuc.Class\tSmall.RI\tSmall.Class\tAny.Class\n")
            else:
                out.write("Position\tMetal.Proba\tMetal.Class\tNuclear.Proba\tNuclear.Class\tSmall.Proba\tSmall.Class"
                          "\tAny.Class\n")
            for idx, p in enumerate(predictions):
                pos = idx + 1

                metal_proba = p[0]
                nuc_proba = p[1]
                small_proba = p[2]

                metal_ri = GeneralInformation.convert_proba_to_ri(metal_proba)
                nuc_ri = GeneralInformation.convert_proba_to_ri(nuc_proba)
                small_ri = GeneralInformation.convert_proba_to_ri(small_proba)

                metal_label = GeneralInformation.get_predicted_label(metal_proba, cutoff)
                nuc_label = GeneralInformation.get_predicted_label(nuc_proba, cutoff)
                small_label = GeneralInformation.get_predicted_label(small_proba, cutoff)

                overall_label = 'nb'
                if metal_label == 'b' or nuc_label == 'b' or small_label == 'b':
                    overall_label = 'b'

                if ri:
                    out.write('{}\t{:.3f}\t{}\t{:.3f}\t{}\t{:.3f}\t{}\t{}\n'.format(pos, metal_ri, metal_label, nuc_ri,
                                                                                    nuc_label, small_ri, small_label,
                                                                                    overall_label))
                else:
                    out.write('{}\t{:.3f}\t{}\t{:.3f}\t{}\t{:.3f}\t{}\t{}\n'.format(pos, metal_proba, metal_label,
                                                                                    nuc_proba, nuc_label, small_proba,
                                                                                    small_label, overall_label))


class GeneralInformation(object):

    @staticmethod
    def get_predicted_label(proba, cutoff):
        if proba >= cutoff:
            return 'b'
        else:
            return 'nb'

    @staticmethod
    def convert_proba_to_ri(proba):
        """Convert probabilitiy to RI ranging from 0 to 9"""

        if proba < 0.5:
            ri = round((0.5 - proba) * 9 / 0.5)
        else:
            ri = round((proba - 0.5) * 9 / 0.5)

        return ri

    @staticmethod
    def seed_worker(worker_id):
        worker_seed = torch.initial_seed() % 2 ** 32
        numpy.random.seed(worker_seed)
        random.seed(worker_seed)

    @staticmethod
    def seed_all(seed):
        if not seed:
            seed = 10

        # print("[ Using Seed : ", seed, " ]")

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.cuda.manual_seed(seed)
        numpy.random.seed(seed)
        random.seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

    @staticmethod
    def remove_padded_positions(pred, target, i):
        indices = (i[i.shape[0] - 1, :] != 0).nonzero()

        pred_i = pred[:, indices].squeeze()
        target_i = target[:, indices].squeeze()

        return pred_i, target_i
