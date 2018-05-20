""" Select a model specified by used input values"""


class SelectModel(object):
    def define_model(self, model):
        cols_to_remove = []
        ranges = []
        words = []

        if model == "mm3_cons_solv":
            cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14,
                              16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39,
                              41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "no"]
        if model == "mm3_cons_solv_more":
            cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14,
                              16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39,
                              41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_cons_solv_no_dist_more":
            cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39,
                              41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_cs_3":
            cols_to_remove = [3, 4, 8, 9, 13, 14,
                              18, 19, 23, 24, 28, 29, 33, 34, 38, 39,
                              43, 44, 48, 49]
            ranges = []
            words = ["no", "no"]
        if model == "mm3_cs_3_more":
            cols_to_remove = [3, 4, 8, 9, 13, 14,
                              18, 19, 23, 24, 28, 29, 33, 34, 38, 39,
                              43, 44, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_cs_5":
            cols_to_remove = []
            ranges = []
            words = ["no", "no"]
        if model == "mm3_cs_5_more":
            cols_to_remove = []
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_3":
            cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39,
                              40, 41, 42, 43, 44, 45, 46, 47, 48, 49]
            ranges = []
            words = ["no", "no"]
        if model == "mm3_3_more":
            cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39,
                              40, 41, 42, 43, 44, 45, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_3_cons":
            cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39,
                              40, 41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "no"]
        if model == "mm3_3_cons_more":
            cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39,
                              40, 41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_3_cons_solv":
            cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39, 41, 42, 43, 44, 46, 47, 48,
                              49]
            ranges = []
            words = ["no", "no"]
        if model == "mm3_3_cons_solv_more":
            cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39, 41, 42, 43, 44, 46, 47, 48,
                              49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_5":
            cols_to_remove = [40, 41, 42, 43, 44, 45, 46, 47, 48, 49]
            ranges = []
            words = ["no", "no"]
        if model == "mm3_5_more":
            cols_to_remove = [40, 41, 42, 43, 44, 45, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_5_cons":
            cols_to_remove = [40, 41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "no"]
        if model == "mm3_5_cons_more":
            cols_to_remove = [40, 41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_5_cons_solv":
            cols_to_remove = [41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "no"]
        if model == "mm3_5_cons_solv_more":
            cols_to_remove = [41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm3_avg_lr":
            cols_to_remove = [40, 41, 42, 43, 44, 45, 46, 47, 48, 49]
            ranges = [[0, 0], [1, 2], [3, 4], [5, 5], [6, 7], [8, 9], [10, 10], [11, 12], [13, 14], [15, 15], [16, 17],
                      [18, 19],
                      [20, 20], [21, 22], [23, 24], [25, 25], [26, 27], [28, 29], [30, 30], [31, 32], [33, 34],
                      [35, 35], [36, 37], [38, 39]]
            words = ["yes", "no"]
        if model == "mm3_avg_lr_more":
            cols_to_remove = [40, 41, 42, 43, 44, 45, 46, 47, 48, 49]
            ranges = [[0, 0], [1, 2], [3, 4], [5, 5], [6, 7], [8, 9], [10, 10], [11, 12], [13, 14], [15, 15], [16, 17],
                      [18, 19],
                      [20, 20], [21, 22], [23, 24], [25, 25], [26, 27], [28, 29], [30, 30], [31, 32], [33, 34],
                      [35, 35], [36, 37], [38, 39]]
            words = ["yes", "yes"]

        if model == "mm_best":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              31, 32, 33, 34, 36, 37, 38, 39, 41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "no"]
        if model == "mm_best_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              31, 32, 33, 34, 36, 37, 38, 39, 41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm_best_3":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14,
                              15, 16, 17, 18, 19, 23, 24, 25, 26, 27, 28, 29,
                              33, 34, 38, 39, 43, 44, 48, 49]
            ranges = []
            words = ["no", "no"]
        if model == "mm_best_3_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14,
                              15, 16, 17, 18, 19, 23, 24, 25, 26, 27, 28, 29,
                              33, 34, 38, 39, 43, 44, 48, 49]
            ranges = []
            words = ["no", "yes"]
        if model == "mm_best_5":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                              15, 16, 17, 18, 19, 25, 26, 27, 28, 29]
            ranges = []
            words = ["no", "no"]
        if model == "mm_best_5_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                              15, 16, 17, 18, 19, 25, 26, 27, 28, 29]
            ranges = []
            words = ["no", "yes"]

        if model == "snap_cons_solv_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              30, 31, 32, 33, 34, 36, 37, 38, 39,
                              41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]

        if model == "snap_cons_solv_3":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              30, 31, 32, 33, 34, 38, 39,
                              43, 44, 48, 49]
            ranges = []
            words = ["no", "no"]

        if model == "snap_cons_solv_3_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              30, 31, 32, 33, 34, 38, 39,
                              43, 44, 48, 49]
            ranges = []
            words = ["no", "yes"]

        if model == "snap_cons_solv_5_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              30, 31, 32, 33, 34]
            ranges = []
            words = ["no", "yes"]

        if model == "evmut_cons_solv_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              31, 32, 33, 34, 35, 36, 37, 38, 39,
                              41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]

        if model == "evmut_cons_solv_3_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              33, 34, 35, 36, 37, 38, 39,
                              43, 44, 48, 49]
            ranges = []
            words = ["no", "yes"]

        if model == "evmut_cons_solv_5_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              35, 36, 37, 38, 39]
            ranges = []
            words = ["no", "yes"]

        if model == "snap_evmut_cons_solv_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              31, 32, 33, 34, 36, 37, 38, 39,
                              41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]

        if model == "snap_evmut_cons_solv_3_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                              33, 34, 38, 39,
                              43, 44, 48, 49]
            ranges = []
            words = ["no", "yes"]

        if model == "snap_evmut_cons_solv_5_more":
            cols_to_remove = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                              15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
            ranges = []
            words = ["no", "yes"]

        return cols_to_remove, ranges, words
