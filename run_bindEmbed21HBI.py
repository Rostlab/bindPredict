import numpy
from pathlib import Path

from config import FileManager, FileSetter
from data_preparation import ProteinInformation, ProteinResults

from annotation_transfer import Inference


def main():

    set_name = 'example'  # TODO set a meaningful set name used for output files
    query_fasta = FileSetter.query_set()
    query_sequences = FileManager.read_fasta(query_fasta)
    query_ids = list(query_sequences.keys())

    result_folder = FileSetter.mmseqs_output()
    Path(result_folder).mkdir(parents=True, exist_ok=True)

    e_val = '1e-3'  # E-value used for bindEmbed21HBI, can be changed
    proteins_inferred = dict()

    # inference set with annotations
    hbi_sequences = FileManager.read_fasta(FileSetter.lookup_fasta())
    labels_hbi = ProteinInformation.get_labels(list(hbi_sequences.keys()), hbi_sequences)

    inference = Inference(query_ids, labels_hbi, query_sequences)
    # run HBI and don't allow self-hits
    inferred_labels, inferred_proteins = inference.infer_binding_annotations_seq(e_val, 'eval', set_name)

    for p in inferred_proteins:
        prot = ProteinResults(p, 3)
        prot.set_predictions(numpy.transpose(inferred_labels[p]))
        proteins_inferred[p] = prot

    print('Number of proteins with hit: {}'.format(len(proteins_inferred.keys())))

    # write predictions
    predictions_folder = '{}/hbi_predictions/'.format(FileSetter.mmseqs_output())
    Path(predictions_folder).mkdir(parents=True, exist_ok=True)
    FileManager.write_predictions(proteins_inferred, predictions_folder, 0.5, True)


main()
