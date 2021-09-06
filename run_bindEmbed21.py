import numpy
from pathlib import Path

from config import FileManager, FileSetter
from data_preparation import ProteinInformation, ProteinResults

from annotation_transfer import Inference
from bindEmbed21DL import BindEmbed21DL


def main():

    set_name = 'example'  # TODO set a meaningful set name used for output files
    query_fasta = FileSetter.query_set()
    query_sequences = FileManager.read_fasta(query_fasta)
    query_ids = list(query_sequences.keys())

    print('Run predictions using homology-based inference (bindEmbed21HBI)')
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

    missing_proteins = []
    for p in query_ids:
        if p in inferred_proteins:
            prot = ProteinResults(p, 3)
            prot.set_predictions(numpy.transpose(inferred_labels[p]))
            proteins_inferred[p] = prot
        else:
            missing_proteins.append(p)

    print('Number of proteins with hit: {}'.format(len(proteins_inferred.keys())))

    print('Run predictions using Machine Learning for remaining proteins (bindEmbed21DL)')
    model_prefix = 'trained_models/checkpoint'
    ri = True  # Whether to write RI or Probabilities

    predicted_proteins = BindEmbed21DL.prediction_pipeline(model_prefix, 0.5, None, missing_proteins, query_sequences,
                                                           ri)

    print("Write predictions")
    proteins = {**inferred_proteins, **predicted_proteins}
    prediction_folder = '{}/predictions'.format(result_folder)
    Path(prediction_folder).mkdir(parents=True, exist_ok=True)
    FileManager.write_predictions(proteins, prediction_folder, 0.5, ri)


main()
