# bindPredict

bindPredict is a sequence-based method to predict binding site residues. 
It uses an Artifical Neural Network (ANN) and was trained on enzymes and DNA-binding proteins 
(i.e. it is most suited to predict binding sites in these type of proteins).
It uses 10 input features calculated for each residue:
- cumulative coupling scores
- cumulative coupling scores filtered by distance
- cumulative coupling scores filtered by solvent accessibility
- clustering coefficients
- clustering coefficients filtered by distance
- clustering coefficients filterd by solvent accessibility
- BLOSUM62-SNAP2 comparison
- BLOSUM62-EVmutation comparison
- relative solvent accessibility
- conservation


## Implementation
bindPredict is written in Python3. The trained model using the 10 input features named above is included in this repository.
bindPredict can be executed via command line runnning "bindPredict.py".
To calculate the different scores used for the prediction, bindPredict needs pre-calculated results from external programs.
As input parameters, the method needs:
- path to a folder containing EVcouplings results 
  (including *.di_scores, *_CouplingScoresCompared_all.csv, *.epistatic_effect, *_frequencies.csv, *_alignment_statistics.csv)
- path to a folder containing SNAP2 results
- path to a folder containing PROFphd results

The method then makes a prediction for every protein present in the SNAP2 results folder. 
If any files are not provided, the corresponding scores are not calculated and set to 0.

## External programs
The pre-calculated results from external programs can be calculated as follows:

- **SNAP2** results: SNAP2 can be run using [predictprotein.org](https://www.predictprotein.org/) [1]. After the job has finished, the results can be downloaded under "Effect of Point Mutations" using the export button.
- **PROFphd** results: PROFphd can also be run using [predictprotein.org](https://www.predictprotein.org/). After the job has finised, the view of the results have to be changed to "text" (upper right corner) and the text part under "Secondary Structure" has to be copied to a text file.
- **EVcouplings** results (\*.di_scores, \*_CouplingScoresCompared_all.csv, \*_frequencies.csv, \*_alignment_statistics.csv): [EVcouplings](http://evfold.org/evfold-web/newmarkec.do) is available as a webserver.
- **EVmutation** results (\*.epistatic_effect): 
[EVmutation](https://marks.hms.harvard.edu/evmutation/index.html) [3] is only available as a Github repository. A detailed description on how to run EVmutation can be found [here](https://htmlpreview.github.io/?https://github.com/debbiemarkslab/EVmutation/blob/master/EVmutation.html)

## References
[1] Yachdav, G., Kloppmann, E., Laszlo, K., Hecht, M., Goldberg, T., Hamp, T., Hönigschmid, P., Schafferhans, A., Roos, M., Bernhofer, M. et al. (2014). [PredictProtein - an open resource for online prediction of protein structural and functional features](https://academic.oup.com/nar/article/42/W1/W337/2435518). Nucleic acids research **42**, pages W337–W343.

[2]EVcouplings

[3] Hopf, T. A., Ingraham, J. B., Poelwijk, F.J., Schärfe, C.P.I., Springer, M., Sander, C., & Marks, D. S. (2017). 
[Mutation effects predicted from sequence co-variation](https://www.nature.com/articles/nbt.3769). Nature Biotechnology **35**, pages 128–135.

