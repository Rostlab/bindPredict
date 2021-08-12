"""
DESCRIPTION:

Development and training of bindPredict to predict binding site residues
"""

import argparse
import sys

from scripts.helper import Helper
from scripts.bindpredict_predictor import Predictor


def main():
    usage_string = 'python bindPredict.py evcouplings_res/ snap_res/ solv_acc/ results/'

    # parse command line options
    parser = argparse.ArgumentParser(description=__doc__, usage=usage_string)
    parser.add_argument('evc_folder', help='Folder with protein folders containing EVcouplings results')
    parser.add_argument('snap_folder', help='Folder with SNAP2 predictions. '
                                            'SNAP-Files need to have the same name in front of the suffix as '
                                            'the EVcouplings folders '
                                            '(e.g. Q9XLZ3/ <-> Q9XLZ3.snap2)')
    parser.add_argument('solv_acc_folder', help='Folder with PROFacc predictions. '
                                                'PROFacc-Files need to have the same name in front of the suffix as '
                                                'the EVcouplings folders '
                                                '(e.g. Q9XLZ3/ <-> Q9XLZ3.reprof)')
    parser.add_argument('output_folder', help='Folder to write binding site predictions to (one file for each protein)')
    parser.add_argument('-m', '--model', help='Input model to use for predictor training (default: mm3_cons_solv_more)',
                        default='mm3_cons_solv_more')
    parser.add_argument('--snap_suffix', help='Suffix of files in given SNAP-folder (default: "*.snap2")',
                        default='*.snap2')
    parser.add_argument('--profacc_suffix', help='Suffix of files in given PROFacc-folder (default: "*.reprof")',
                        default='*.reprof')
    args = parser.parse_args()
    helper = Helper()
    if helper.folder_existence_check(args.evc_folder):  # check if evc folder exists and is reachable
        if helper.folder_existence_check(args.snap_folder):  # check if snap folder exists and is reachable
            if helper.folder_existence_check(
                    args.solv_acc_folder):  # check if solvent accessibility folder exists and is reachable
                predictor = Predictor()
                predictor.run_prediction(args)


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)


if __name__ == "__main__":
    main()
