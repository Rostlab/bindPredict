import os

class Helper(object):
    def folder_existence_check(self, folder_name=None):
        """
        Check if folder exists and is reachable, if not exit the program.
        :param folder_name: Folder to check, if it exists (default: None)
        :type folder_name: String
        :return: True if folder exists and is reachable, False otherwise
        """
        folder_is_there = False
        if os.path.exists(folder_name):
            folder_is_there = True
        else:
            error('Folder {fl} does not exist or is not reachable - exit!'.format(fl=folder_name))
            folder_is_there = False
            exit(404)

        return folder_is_there


    def read_id_lists(self, ids_folder, ids_file):
        """
        Read the id lists from a given folder
        :param ids_folder: Folder were ids files are located
        :return: ids of split 1-5
        """

        split1_ids = []
        split2_ids = []
        split3_ids = []
        split4_ids = []
        split5_ids = []
        
        with open(os.path.join(ids_folder, (ids_file + "1.txt"))) as f:
            for line in f:
                split1_ids.append(line.strip())

        with open(os.path.join(ids_folder, (ids_file + "2.txt"))) as f:
            for line in f:
                split2_ids.append(line.strip())

        with open(os.path.join(ids_folder, (ids_file + "3.txt"))) as f:
            for line in f:
                split3_ids.append(line.strip())

        with open(os.path.join(ids_folder, (ids_file + "4.txt"))) as f:
            for line in f:
                split4_ids.append(line.strip())

        with open(os.path.join(ids_folder, (ids_file + "5.txt"))) as f:
            for line in f:
                split5_ids.append(line.strip())

        return split1_ids, split2_ids, split3_ids, split4_ids, split5_ids


    def read_more_data(self, ids_folder, run):
        """
        Read ids for more data for the training set of the specified run
        :param ids_folder: Folder were id files are located
        :return: ids for more data
        """
        
        ids = []
        with open(os.path.join(ids_folder, ("ids_run" + str(run) + "_more.txt"))) as f:
            for line in f:
                ids.append(line.strip())

        return ids

