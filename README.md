# bindPredictDL21

bindPredictDL21 is a method to predict whether a residue in a protein is binding to metal ions, nucleic acids (DNA or RNA), or small molecules. Towards this end, bindPredictDL21 combines homology-based inference and Machine Learning. Homology-based inference is executed using MMseqs2 [1]. For the Machine Learning method, bindPredictDL21 uses ProtT5 embeddings [2] as input to a 2-layer CNN. Since bindPredictDL21 is based on single sequences, it can easily be applied to any protein sequence.

## How to use bindPredictDL21

## Data

The data set used for training and testing was extracted from BioLip [3]. The UniProt identifiers for the 5 splits used during cross-validation and the test set as well as the corresponding FASTA sequences and used binding annotations are made available in the `data` folder.

The trained models are available in the `trained_models` folder.

ProtT5 embeddings can be generated using [the bio_embeddings pipeline](https://github.com/sacdallago/bio_embeddings) [4].

## Requirements

bindPredictDL21 is written in Python3. In order to execute bindPredictDL21, Python3 has to be installed locally. Additionally, the following Python packages have to be installed:
- numpy
- scikit-learn
- torch
- pandas

To be able to run homology-based inference, MMseqs2 has to be locally installed. Otherwise, it is also possible to only run the Machine Learning part of bindPredictDL21.

## References
[1] Steinegger M, Söding J (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35.

[2] Elnaggar A, Heinzinger M, Dallago C, Rihawi G, Wang Y, Jones L, Gibbs T, Feher T, Angerer C, Bhowmik D, Rost B (2021). ProtTrans: towards cracking the language of life's code through self-supervised deep learning and high performance computing. bioRxiv.

[3] Yang J, Roy A, Zhang Y (2013). BioLip: a semi-manually curated database for biologically relevant ligand-protein interactions. Nucleic Acids Research, 41.

[4] Dallago C, Schütze K, Heinzinger M, Olenyi T, Littmann M, Lu AX, Yang KK, Min S, Yoon S, Morton JT, & Rost B (2021). Learned embeddings from deep learning to visualize and predict protein sets. Current Protocols, 1, e113. doi: 10.1002/cpz1.113


## bindPredictML17
If you are interested in the predecessor of bindPredictDL21, bindPredictML17, you can find all relevant information in the subfolder `bindPredictML17`.
