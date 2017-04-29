import random
import train_model as tm

# generate id lists for train, cross-train and test splits

ids_dna = []
ids_enzyme = []
ids_more_data = []

with open("data/ids_dna.txt") as f:
    for row in f:
        ids_dna.append(row.strip())

with open("data/ids_enzyme.txt") as f:
    for row in f:
        ids_enzyme.append(row.strip())

with open("data/ids_more_data.txt") as f:
    for row in f:
        ids_more_data.append(row.strip())

random_indices_dna = random.sample(range(0, len(ids_dna)), 5)
random_indices_enzyme = random.sample(range(0, len(ids_enzyme)), 36)

test_dna = [ids_dna[x] for x in random_indices_dna]
test_enzyme = [ids_enzyme[x] for x in random_indices_enzyme]

ids_dna = [i for j, i in enumerate(ids_dna) if j not in random_indices_dna]
ids_enzyme = [i for j, i in enumerate(ids_enzyme) if j not in random_indices_enzyme]

ids = ids_dna + ids_enzyme
random_indices = random.sample(range(0, len(ids)), len(ids))

split1_indices = random_indices[0:74]
split2_indices = random_indices[75:149]
split3_indices = random_indices[150:224]
split4_indices = random_indices[225:298]
split5_indices = random_indices[299:372]

test_ids = test_dna + test_enzyme
split1_ids = [ids[x] for x in split1_indices]
split2_ids = [ids[x] for x in split2_indices]
split3_ids = [ids[x] for x in split3_indices]
split4_ids = [ids[x] for x in split4_indices]
split5_ids = [ids[x] for x in split5_indices]

# define models, calculate accuracy and prec-cov-curve

model = "mm3"
print(model)
cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14,
                  16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39, 40, 41]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "yes"])

model = "mm3_more"
print(model)
cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14,
                  16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39, 40, 41]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "yes"])

model = "mm3_cons"
print(model)
cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14,
                  16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39, 40]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "no"])

model = "mm3_cons_more"
print(model)
cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14,
                  16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39, 40]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "no"])

model = "mm3_cons_solv"
print(model)
cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14,
                  16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "no"])

model = "mm3_cons_solv_more"
print(model)
cols_to_remove = [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14,
                  16, 17, 18, 19, 21, 22, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 36, 37, 38, 39]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "no"])

model = "mm3_3"
print(model)
cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39, 40, 41]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "yes"])

model = "mm3_3_more"
print(model)
cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39, 40, 41]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "yes"])

model = "mm3_3_cons"
print(model)
cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39, 40]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "no"])

model = "mm3_3_cons_more"
print(model)
cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39, 40]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "no"])

model = "mm3_3_cons_solv"
print(model)
cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "no"])

model = "mm3_3_cons_solv_more"
print(model)
cols_to_remove = [3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "no"])

model = "mm3_5"
print(model)
cols_to_remove = [40, 41]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "yes"])

model = "mm3_5_more"
print(model)
cols_to_remove = [40, 41]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "yes"])

model = "mm3_5_cons"
print(model)
cols_to_remove = [40]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "no"])

model = "mm3_5_cons_more"
print(model)
cols_to_remove = [40]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "no"])

model = "mm3_5_cons_solv"
print(model)
cols_to_remove = []
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "no"])

model = "mm3_5_cons_solv_more"
print(model)
cols_to_remove = []
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "no"])

model = "mm3_avg_lr"
print(model)
cols_to_remove = [40, 41]
ranges = [[0,0], [1,2], [3,4], [5,5], [6,7], [8,9], [10,10], [11,12], [13,14], [15,15], [16,17], [18,19],
          [20,20], [21,22], [23,24], [25,25], [26,27], [28,29], [30,30], [31,32], [33,34], [35,35], [36,37], [38,39]]

tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["yes", "no", "no"])

model = "mm3_avg_lr_more"
print(model)
cols_to_remove = [40, 41]
ranges = [[0,0], [1,2], [3,4], [5,5], [6,7], [8,9], [10,10], [11,12], [13,14], [15,15], [16,17], [18,19],
          [20,20], [21,22], [23,24], [25,25], [26,27], [28,29], [30,30], [31,32], [33,34], [35,35], [36,37], [38,39]]

tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["yes", "yes", "no"])

model = "mm3_avg"
print(model)
cols_to_remove = [40, 41]
ranges = [[0,0], [1,4], [5,5], [6,9], [10,10], [11,14], [15,15], [16,19],
          [20,20], [21,24], [25,25], [26,29], [30,30], [31,34], [35,35], [36,39]]

tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["yes", "no", "no"])

model = "mm3_avg_more"
print(model)
cols_to_remove = [40, 41]
ranges = [[0,0], [1,4], [5,5], [6,9], [10,10], [11,14], [15,15], [16,19],
          [20,20], [21,24], [25,25], [26,29], [30,30], [31,34], [35,35], [36,39]]

tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["yes", "yes", "no"])

model = "snap_center"
print(model)
cols_to_remove = [40, 41]
ranges = [[0,4], [5,9], [10,14], [15,19],
          [20,24], [25,29], [30,34], [35,35], [36,39]]
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["yes", "no", "no"])

model = "snap_center_more"
print(model)
cols_to_remove = [40, 41]
ranges = [[0,4], [5,9], [10,14], [15,19],
          [20,24], [25,29], [30,34], [35,35], [36,39]]
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["yes", "yes", "no"])

model = "snap_evc_center"
print(model)
cols_to_remove = [40, 41]
ranges = [[0,4], [5,9], [10,14], [15,19],
          [20,24], [25,29], [30,30], [31,34], [35,35], [36,39]]
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["yes", "no", "no"])

model = "snap_evc_center_more"
print(model)
cols_to_remove = [40, 41]
ranges = [[0,4], [5,9], [10,14], [15,19],
          [20,24], [25,29], [30,30], [31,34], [35,35], [36,39]]
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["yes", "yes", "no"])

model = "mm3_5_wt_cum"
print(model)
cols_to_remove = [0, 1, 2, 3, 4, 40, 41]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "no", "no"])

model = "mm3_5_wt_cum_more"
print(model)
cols_to_remove = [0, 1, 2, 3, 4, 40, 41]
ranges = []
tm.run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
             words=["no", "yes", "no"])
