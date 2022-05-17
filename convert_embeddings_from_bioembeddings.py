"""
Script to convert embeddings generated with the bio_embeddings pipeline in the expected format
"""

import h5py
import numpy as np
import sys


def main():
    embeddings_in = sys.argv[1]
    embeddings_out = sys.argv[2]

    with h5py.File(embeddings_in, 'r') as f_in:
        with h5py.File(embeddings_out, 'w') as f_out:
            for key,embedding in f_in.items():
                original_id = embedding.attrs['original_id']
                embedding = np.array(embedding)

                f_out.create_dataset(original_id, data=embedding)


main()
