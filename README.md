# bindEmbed21

bindEmbed21 is a method to predict whether a residue in a protein is binding to metal ions, nucleic acids (DNA or RNA), or small molecules. Towards this end, bindEmbed21 combines homology-based inference and Machine Learning. Homology-based inference is executed using MMseqs2 [1]. For the Machine Learning method, bindEmbed21DL uses ProtT5 embeddings [2] as input to a 2-layer CNN. Since bindEmbed21 is based on single sequences, it can easily be applied to any protein sequence.

## Usage

`run_bindEmbed21DL.py` shows an example how to generate binding residue predictions using the Machine Learning part of bindEmbed21 (bindEmbed21DL)

`run_bindEmbed21HBI.py` shows an example how to generate bidning residue predictions using the homology-inference part of bindEmbed21 (bindEmbed21HBI)

`run_bindEmbed21.py` combines ML and HBI into the final method bindEmbed21

`develop_bindEmbed21DL.py` provides the code to reproduce the bindEmbed21DL development (hyperparameter optimization, training, performance assessment on the test set).

All needed files and paths can be set in `config.py` (marked as TODOs).

## Data

### Development Set

The data set used for training and testing was extracted from BioLip [3]. The UniProt identifiers for the 5 splits used during cross-validation (DevSet1014), the test set (TestSet300), and the independent set of proteins added to BioLip after November 2019 (TestSetNew46) as well as the corresponding FASTA sequences and used binding annotations are made available in the `data` folder.

The trained models are available in the `trained_models` folder.

ProtT5 embeddings can be generated using [the bio_embeddings pipeline](https://github.com/sacdallago/bio_embeddings) [4]. To use them with `bindEmbed21`, they need to be converted to use the correct keys. A script for the conversion can be found in the folder `utils`.

## Sets for homology-based inference
For the homology-based inference (bindEmbed21HBI), query proteins will be aligned against big80 to generate profiles. Those profiles are then searched against a lookup set of proteins with known binding residues. The pre-computed MMseqs2 database files and the FASTA file for the lookup database can be downloaded here:

* Pre-computed big80 DB: [ftp://rostlab.org/bindEmbed21/profile_db.tar.gz](ftp://rostlab.org/bindEmbed21/profile_db.tar.gz)
* Pre-computed lookup DB: [ftp://rostlab.org/bindEmbed21/lookup_db.tar.gz](ftp://rostlab.org/bindEmbed21/lookup_db.tar.gz)
* FASTA for lookup DB: [ftp://rostlab.org/bindEmbed21/lookup.fasta](ftp://rostlab.org/bindEmbed21/lookup.fasta)

### Human proteome predictions

We applied bindEmbed21DL as well as homology-based inference to the entire human proteome. While annotations were only available for 15% of the human proteins, homology-based inference allowed transferring annotations for 48% (9,694) and bindEmbed21DL provided binding predictions for 92% (18,663) of the human proteome. Both predictions are available in the folder `human_proteome`. For predictions made using homology-based inference, values of -1.0 refer to position which were not inferred, and therefore, were considered non-binding.

## Requirements

bindEmbed21 is written in Python3. In order to execute bindEmbed21, Python3 has to be installed locally. Additionally, the following Python packages have to be installed:
- numpy
- scikit-learn
- torch
- pandas
- h5py

To be able to run homology-based inference, MMseqs2 has to be locally installed. Otherwise, it is also possible to only run the Machine Learning part of bindEmbed21 (bindEmbed21DL).

## Cite

In case, you are using this method and find it helpful, we would appreciate if you could cite the following publication:

Littmann M, Heinzinger M, Dallago C, Weissenow K, Rost B. Protein embeddings and deep learning predict binding residues for various ligand classes. *Sci Rep* **11**, 23916 (2021). https://doi.org/10.1038/s41598-021-03431-4


## References
[1] Steinegger M, Söding J (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35.

[2] Elnaggar A, Heinzinger M, Dallago C, Rihawi G, Wang Y, Jones L, Gibbs T, Feher T, Angerer C, Bhowmik D, Rost B (2021). ProtTrans: towards cracking the language of life's code through self-supervised deep learning and high performance computing. bioRxiv.

[3] Yang J, Roy A, Zhang Y (2013). BioLip: a semi-manually curated database for biologically relevant ligand-protein interactions. Nucleic Acids Research, 41.

[4] Dallago C, Schütze K, Heinzinger M, Olenyi T, Littmann M, Lu AX, Yang KK, Min S, Yoon S, Morton JT, & Rost B (2021). Learned embeddings from deep learning to visualize and predict protein sets. Current Protocols, 1, e113. doi: 10.1002/cpz1.113


## bindPredictML17
If you are interested in the predecessor of bindEmbed21, bindPredictML17, you can find all relevant information in the subfolder `bindPredictML17`.
