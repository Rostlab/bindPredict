# bindPredictDL21

bindPredictDL21 is a method to predict whether a residue in a protein is binding to metal ions, nucleic acids (DNA or RNA), or small molecules. Towards this end, bindPredictDL21 uses ProtT5 embeddings [1] as input to a 2-layer CNN. Since bindPredictDL21 is based on single sequences, it can easily be applied to any protein sequence.

## Data

The data set used for training and testing was extracted from BioLip [2]. The UniProt identifiers for the 5 splits used during cross-validation and the test set as well as the corresponding FASTA sequences and used binding annotations are made available in the `data` folder.

ProtT5 embeddings can be generated using [the bio_embeddings pipeline](https://github.com/sacdallago/bio_embeddings) [3].

## References

[1] Elnaggar A, Heinzinger M, Dallago C, Rihawi G, Wang Y, Jones L, Gibbs T, Feher T, Angerer C, Bhowmik D, Rost B (2021). ProtTrans: towards cracking the language of life's code through self-supervised deep learning and high performance computing. bioRxiv.

[2] Yang J, Roy A, Zhang Y (2013). BioLip: a semi-manually curated database for biologically relevant ligand-protein interactions. Nucleic Acids Research, 41.

[3] Dallago C, Schütze K, Heinzinger M, Olenyi T, Littmann M, Lu AX, Yang KK, Min S, Yoon S, Morton JT, & Rost B (2021). Learned embeddings from deep learning to visualize and predict protein sets. Current Protocols, 1, e113. doi: 10.1002/cpz1.113


## bindPredictML17
If you are interested in the predecessor of bindPredictDL21, bindPredictML17, you can find all relevant information in the subfolder `bindPredictML17`.
