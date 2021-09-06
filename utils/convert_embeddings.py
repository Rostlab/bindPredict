import h5py
import numpy as np


def main():

    embeddings_in = ''  # TODO give path to .h5-file of embeddings generated with bio_embeddings pipeline
    embeddings_out = ''  # TODO set output directory

    embeddings = dict()
    with h5py.File(embeddings_in, 'r') as f:
        for key, embedding in f.items():
            original_id = embedding.attrs['original_id']
            embeddings[original_id] = np.array(embedding)

    for k in embeddings.keys():
        embedding = embeddings[k]

        out_file = '{}{}.npy'.format(embeddings_out, k)
        np.save(out_file, embedding)


main()
