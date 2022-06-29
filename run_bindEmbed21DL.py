from config import FileSetter, FileManager
from bindEmbed21DL import BindEmbed21DL

from pathlib import Path


def main():

    prediction_folder = FileSetter.predictions_folder()
    Path(prediction_folder).mkdir(parents=True, exist_ok=True)

    model_prefix = 'trained_models/checkpoint'

    query_fasta = FileSetter.query_set()
    query_sequences = FileManager.read_fasta(query_fasta)
    query_ids = list(query_sequences.keys())

    ri = True  # Whether to write RI or Probabilities

    BindEmbed21DL.prediction_pipeline(model_prefix, 0.5, prediction_folder, query_ids, query_fasta, ri)


main()
