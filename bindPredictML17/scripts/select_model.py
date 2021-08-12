class SelectModel(object):
    @staticmethod
    def define_model(model):
        cols_to_remove = []
        ranges = []
        words = []

        if model == "mm3_cons_solv_more":
            cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14,
                              16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39,
                              41, 42, 43, 44, 46, 47, 48, 49]
            ranges = []
            words = ["no", "yes"]
        else:
            print("Unknown model specified")

        return cols_to_remove, ranges, words
