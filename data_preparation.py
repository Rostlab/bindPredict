import torch
import numpy as np
import sys
from collections import defaultdict

from config import FileManager, FileSetter
from assess_performance import PerformanceAssessment


class MyDataset(torch.utils.data.Dataset):

    """Dataset for bindEmbed21DL"""

    def __init__(self, samples, seqs, labels, max_length):
        self.samples = samples
        self.seqs = seqs

        self.labels = labels
        self.max_length = max_length

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, item):
        prot_id = self.samples[item]
        prot_length = len(self.seqs[prot_id])

        embedding = np.load('{}/{}.npy'.format(FileSetter.t5_dir(), prot_id))

        # pad all inputs to the maximum length & add another feature to encode whether the element is a position
        # in the sequence or padded
        features = np.zeros((1025, self.max_length), dtype=np.float32)
        features[:1024, :prot_length] = np.transpose(embedding)  # set feature maps to embedding values
        features[1024, :prot_length] = 1  # set last element to 1 because positions are not padded

        target = np.zeros((3, self.max_length), dtype=np.float32)
        target[:3, :prot_length] = np.transpose(self.labels[prot_id])
        loss_mask = np.zeros((3, self.max_length), dtype=np.float32)
        loss_mask[:3, :prot_length] = 1 * prot_length

        if self.protein_prediction:
            return features, target, loss_mask, prot_id
        else:
            return features, target, loss_mask


class ProteinInformation(object):

    @staticmethod
    def get_data(ids):
        """
        Get sequences, labels, and maximum length for a set of ids
        :param ids:
        :return:
        """

        sequences = FileManager.read_fasta(FileSetter.fasta_file())
        max_length = ProteinInformation.determine_max_length(sequences, ids)
        labels = ProteinInformation.get_labels(ids, sequences)

        return sequences, max_length, labels

    @staticmethod
    def get_data_predictions(ids, fasta_file):
        """
        Generate dummy labels for test proteins without annotations to allow re-use of general DataLoader
        :param ids:
        :param fasta_file:
        :return: sequences, max. length, dummy labels
        """
        sequences = FileManager.read_fasta(fasta_file)

        max_length = ProteinInformation.determine_max_length(sequences, ids)
        labels = dict()
        for i in ids:
            prot_length = len(sequences[i])
            binding_tensor = np.zeros([prot_length, 3], dtype=np.float32)
            labels[i] = binding_tensor

        return sequences, max_length, labels

    @staticmethod
    def determine_max_length(sequences, ids):
        """Get maximum length in set of sequences"""
        max_len = 0
        for i in ids:
            if len(sequences[i]) > max_len:
                max_len = len(sequences[i])

        return max_len

    @staticmethod
    def get_labels(ids, sequences, file_prefix=None):
        """
        Read binding residues for metal, nucleic acids, and small molecule binding
        :param ids:
        :param sequences:
        :param file_prefix: If None, files set in FileSetter will be used
        :return:
        """
        labels = dict()

        if file_prefix is None:
            metal_residues = FileManager.read_binding_residues(FileSetter.binding_residues_by_ligand('metal'))
            nuclear_residues = FileManager.read_binding_residues(FileSetter.binding_residues_by_ligand('nuclear'))
            small_residues = FileManager.read_binding_residues(FileSetter.binding_residues_by_ligand('small'))
        else:
            metal_residues = FileManager.read_binding_residues('{}_metal.txt'.format(file_prefix))
            nuclear_residues = FileManager.read_binding_residues('{}_nuclear.txt'.format(file_prefix))
            small_residues = FileManager.read_binding_residues('{}_small.txt'.format(file_prefix))

        for prot_id in ids:
            prot_length = len(sequences[prot_id])
            binding_tensor = np.zeros([prot_length, 3], dtype=np.float32)

            metal_res = nuc_res = small_res = []

            if prot_id in metal_residues.keys():
                metal_res = metal_residues[prot_id]
            if prot_id in nuclear_residues.keys():
                nuc_res = nuclear_residues[prot_id]
            if prot_id in small_residues.keys():
                small_res = small_residues[prot_id]

            metal_residues_0_ind = ProteinInformation._get_zero_based_residues(metal_res)
            nuc_residues_0_ind = ProteinInformation._get_zero_based_residues(nuc_res)
            small_residues_0_ind = ProteinInformation._get_zero_based_residues(small_res)

            binding_tensor[metal_residues_0_ind, 0] = 1
            binding_tensor[nuc_residues_0_ind, 1] = 1
            binding_tensor[small_residues_0_ind, 2] = 1

            labels[prot_id] = binding_tensor

        return labels

    @staticmethod
    def _get_zero_based_residues(residues):
        residues_0_ind = []
        for r in residues:
            residues_0_ind.append(int(r) - 1)

        return residues_0_ind


class ProteinResults(object):

    def __init__(self, name, bind_cutoff=0.5):
        self.name = name
        self.labels = np.array([])
        self.predictions = np.array([])

        self.bind_cutoff = bind_cutoff
        # cutoff to define if label is binding/non-binding; default: 0: non-binding, 1:binding

    def set_labels(self, labels):
        self.labels = np.array(labels)

    def set_predictions(self, predictions):
        self.predictions = np.around(np.array(np.transpose(predictions)), 3)

    def add_predictions(self, predictions):
        self.predictions = np.add(self.predictions, np.around(predictions, 3))

    def normalize_predictions(self, norm_factor):
        self.predictions = np.around(self.predictions / norm_factor, 3)

    def calc_num_predictions(self, cutoff):
        num_predictions = np.count_nonzero(self.predictions >= cutoff, axis=0)

        return num_predictions[0], num_predictions[1], num_predictions[2]

    def get_bound_ligand(self, cutoff):
        num_labels = np.count_nonzero(self.labels >= cutoff, axis=0)

        metal = nuclear = small = False
        if num_labels[0] > 0:
            metal = True
        if num_labels[1] > 0:
            nuclear = True
        if num_labels[2] > 0:
            small = True

        return metal, nuclear, small

    def get_predictions_ligand(self, ligand):

        if ligand == 'metal':
            return self.predictions[:, 0]
        elif ligand == 'nucleic':
            return self.predictions[:, 1]
        elif ligand == 'small':
            return self.predictions[:, 2]
        elif ligand == 'overall':
            return np.amax(self.predictions, axis=1)
        else:
            sys.exit('{} is not a valid ligand type'.format(ligand))

    def calc_performance_measurements(self, cutoff):
        performance = self.calc_performance_given_labels(cutoff, self.labels)

        return performance

    def calc_performance_given_labels(self, cutoff, ligand_labels):
        """Calculate performance values for this protein"""
        performance = defaultdict(dict)
        num_ligands = np.shape(ligand_labels)[1]

        if num_ligands > 1:  # ligand-type assessment
            # calc per-ligand assessment for multi-label prediction
            for i in range(0, num_ligands):
                tp = fp = tn = fn = 0
                cross_pred = [0, 0, 0, 0]
                for idx, lig in enumerate(ligand_labels):
                    if self.predictions[idx, i] >= cutoff:  # predicted as binding to this ligand
                        cross_prediction = False
                        true_prediction = False

                        for j in range(0, num_ligands):
                            if i == j:  # same as predicted ligand
                                if lig[j] >= self.bind_cutoff:  # also annotated to this ligand
                                    tp += 1
                                    cross_pred[i] += 1
                                    true_prediction = True
                                else:
                                    fp += 1
                            else:
                                if lig[j] >= self.bind_cutoff and not true_prediction:
                                    cross_pred[j] += 1
                                    cross_prediction = True

                        if not true_prediction and not cross_prediction:
                            # residues is not annotated to bind any of the ligands
                            cross_pred[3] += 1
                    else:
                        if lig[i] >= cutoff:
                            fn += 1
                        else:
                            tn += 1

                if i == 0:
                    ligand = 'metal'
                elif i == 1:
                    ligand = 'nucleic'
                else:
                    ligand = 'small'

                bound = False
                if (tp + fn) > 0:
                    bound = True
                acc, prec, recall, f1, mcc = PerformanceAssessment.calc_performance_measurements(tp, fp, tn, fn)
                # calculate performance measurements for negatives
                _, neg_p, neg_r, neg_f1, _ = PerformanceAssessment.calc_performance_measurements(tn, fn, tp, fp)

                performance[ligand] = {'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn, 'acc': acc, 'prec': prec,
                                       'recall': recall, 'f1': f1, 'neg_prec': neg_p, 'neg_recall': neg_r,
                                       'neg_f1': neg_f1, 'mcc': mcc, 'bound': bound,
                                       'cross_prediction': cross_pred}

        # get overall performance
        reduced_labels = np.sum(ligand_labels > cutoff, axis=1)
        if len(self.predictions.shape) == 1:
            reduced_predictions = (self.predictions >= cutoff)
        else:
            reduced_predictions = np.sum(self.predictions >= cutoff, axis=1)

        tp = np.sum(np.logical_and(reduced_labels > 0, reduced_predictions > 0))
        fp = np.sum(np.logical_and(reduced_labels == 0, reduced_predictions > 0))
        tn = np.sum(np.logical_and(reduced_labels == 0, reduced_predictions == 0))
        fn = np.sum(np.logical_and(reduced_labels > 0, reduced_predictions == 0))

        acc, prec, recall, f1, mcc = PerformanceAssessment.calc_performance_measurements(tp, fp, tn, fn)
        _, neg_p, neg_r, neg_f1, _ = PerformanceAssessment.calc_performance_measurements(tn, fn, tp, fp)
        performance['overall'] = {'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn, 'acc': acc, 'prec': prec, 'recall': recall,
                                  'f1': f1, 'neg_prec': neg_p, 'neg_recall': neg_r, 'neg_f1': neg_f1, 'mcc': mcc,
                                  'bound': True, 'cross_prediction': [0, 0, 0, 0]}

        return performance
