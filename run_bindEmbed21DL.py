from config import FileSetter, FileManager
from bindEmbed21DL import BindEmbed21DL


def main():

    prediction_folder = FileSetter.predictions_folder()

    model_prefix = 'trained_models/checkpoint'
    predictor = BindEmbed21DL('t5')

    query_fasta = FileSetter.query_set()
    query_sequences = FileManager.read_fasta(query_fasta)
    query_ids = list(query_sequences.keys())

    predictor.prediction_pipeline(model_prefix, 0.5, prediction_folder, query_ids, query_sequences)


main()
