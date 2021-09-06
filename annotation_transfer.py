from homology_based_inference import MMseqsPredictor
from config import FileManager

import numpy
import random


class Inference(object):

    def __init__(self, prot_ids, labels, sequences, rand_pred=False):
        self.ids = prot_ids
        self.labels = labels
        self.sequences = sequences

        self.rand_pred = rand_pred
        self.rnd_probas = {'metal': 0.147, 'nuclear': 0.030, 'small': 0.090}

    def infer_binding_annotations_seq(self, e_val, criterion, set_name, exclude_self_hits=True):

        # get hits from MMseqs2 alignments
        mmseqs_pred = MMseqsPredictor(exclude_self_hits)
        hits, mmseqs_file = mmseqs_pred.get_mmseqs_hits(self.ids, e_val, criterion, set_name)

        inferred_labels, proteins_with_hit = self._infer_binding_local_alignment(hits, mmseqs_file)

        return inferred_labels, proteins_with_hit

    def _infer_binding_local_alignment(self, hits, mmseqs_file):
        """
        Get binding annotations from local alignment
        :param hits: Found hits
        :param mmseqs_file: File with local alignments from MMseqs2
        :return:
        """

        inferred_labels = dict()
        proteins_with_hit = dict()

        mmseqs = FileManager.read_mmseqs_alignments(mmseqs_file)

        for i in hits.keys():
            hit = list(hits[i].keys())[0]

            alignment = mmseqs[i][hit]
            indices1, indices2 = Inference._get_index_mapping_mmseqs(alignment)

            binding_annotation = self.labels[hit]
            # set all position to -1; replace inferred positions with 0 (non-binding) or 1 (binding)
            inferred_annotation = numpy.zeros([len(self.sequences[i]), 3], dtype=numpy.float32) - 1

            for idx, pos2 in enumerate(indices2):
                pos1 = indices1[idx] - 1

                if pos1 >= 0 and pos2 >= 1:  # both positions are aligned
                    anno = binding_annotation[pos2 - 1]
                    inferred_annotation[pos1] = anno

            if 1 in inferred_annotation:
                # any binding annotation in the aligned region
                inferred_labels[i] = inferred_annotation
                ligands = {'metal': 0, 'nuclear': 0, 'small': 0}

                if 1 in inferred_annotation[:, 0]:
                    ligands['metal'] = 1
                if 1 in inferred_annotation[:, 1]:
                    ligands['nuclear'] = 1
                if 1 in inferred_annotation[:, 2]:
                    ligands['small'] = 1

                proteins_with_hit[i] = ligands

        for i in self.ids:
            if i not in inferred_labels.keys():  # no hit generated
                inferred_annotation = numpy.zeros([len(self.sequences[i]), 3], dtype=numpy.float32) - 1

                # generate random prediction for proteins without a hit
                if self.rand_pred:

                    metal_indices = random.choices([0,  1],
                                                   weights=[1-self.rnd_probas['metal'], self.rnd_probas['metal']],
                                                   k=len(self.sequences[i]))
                    nuclear_indices = random.choices([0, 1],
                                                     weights=[1-self.rnd_probas['nuclear'], self.rnd_probas['nuclear']],
                                                     k=len(self.sequences[i]))
                    small_indices = random.choices([0, 1],
                                                   weights=[1-self.rnd_probas['small'], self.rnd_probas['small']],
                                                   k=len(self.sequences[i]))

                    inferred_annotation[:, 0] = metal_indices
                    inferred_annotation[:, 1] = nuclear_indices
                    inferred_annotation[:, 2] = small_indices

                inferred_labels[i] = inferred_annotation

        return inferred_labels, proteins_with_hit

    @staticmethod
    def _get_index_mapping_mmseqs(alignment):
        """
        Get indices of alignment position in actual sequence
        :param alignment: Local alignment of 2 sequences
        :return: Indices for 2 sequences
        """

        start1 = alignment['qstart']
        start2 = alignment['tstart']

        seq1 = alignment['qaln']
        seq2 = alignment['taln']

        indices1 = Inference._get_indices_seq(start1, seq1)
        indices2 = Inference._get_indices_seq(start2, seq2)

        return indices1, indices2

    @staticmethod
    def _get_indices_seq(start, seq):
        """
        Get sequence indices of aligned sequence
        :param start: Start index of aligned sequence
        :param seq: Aligned sequence (incl. gaps)
        :return: Sequence indices for position in aligned sequence
        """

        indices = []

        for i in range(0, len(seq)):
            if seq[i] == '-':  # position is gap
                indices.append(0)
            else:
                indices.append(start)
                start += 1

        return indices
