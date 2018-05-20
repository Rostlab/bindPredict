""" Class specifying file endings and locations.
Before usage, values should be adjusted to local settings"""


class FileManager(object):

    @property
    def ec_suffix(self):
        return ".di_scores"

    @property
    def dist_suffix(self):
        return "_CouplingScoresCompared_all.csv"

    @property
    def evmut_suffix(self):
        return ".epistatic_effect"

    @property
    def cons_suffix(self):
        return "_frequencies.csv"

    @property
    def align_suffix(self):
        return "_alignment_statistics.csv"

    @property
    def ids_file(self):
        file_name = "ids_split"
        return file_name

    @property
    def ids_folder(self):
        return "data/id_lists/"

    @property
    def binding_file(self):
        return "data/binding_residues_evc.txt"

    @property
    def blosum_file(self):
        return "data/blosum_62.txt"
