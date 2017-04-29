import functions as func
from imblearn.over_sampling import SMOTE
from sklearn.neural_network import MLPClassifier
from sklearn import metrics
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt


def run_model(model, split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
              words):

    x_train, y_train, train_length, x_cross, y_cross, cross_length = import_model(split1_ids, split2_ids, split3_ids,
                                                                                  split4_ids, split5_ids,
                                                                                  ids_more_data, cols_to_remove,
                                                                                  ranges, keywords=words[0:1])
    classifier = train_model(x_train, y_train)
    mean_line, sd_line = model_accuracy(classifier, x_cross, y_cross, cross_length)
    prec_cov = prec_cov_curve(classifier, x_cross, y_cross, cross_length)

    mean_line = model + mean_line
    sd_line = model + sd_line

    with open("data/performance/performance_mean_new.txt", "a") as f:
        f.write(mean_line)
    with open("data/performance/performance_sd_new.txt", "a") as f:
        f.write(sd_line)
    with open("data/prec_cov/prec_cov_" + model + ".txt", "w") as f:
        for x in prec_cov[0]:
            f.write(x + "\t")
        f.write("\n")
        for x in prec_cov[1]:
            f.write(x + "\t")
        f.write("\n")

    if words[2] == "yes":
        test_ovetraining(model,x_train, y_train, x_cross, y_cross)


def import_model(split1_ids, split2_ids, split3_ids, split4_ids, split5_ids, ids_more_data, cols_to_remove, ranges,
                 keywords):
    r = keywords[0]
    more = keywords[1]

    x_split1, y_split1, length_split1 = func.import_data(split1_ids, cols_to_remove)
    x_split2, y_split2, length_split2 = func.import_data(split2_ids, cols_to_remove)
    x_split3, y_split3, length_split3 = func.import_data(split3_ids, cols_to_remove)
    x_split4, y_split4, length_split4 = func.import_data(split4_ids, cols_to_remove)
    x_split5, y_split5, length_split5 = func.import_data(split5_ids, cols_to_remove)
    x_more, y_more, length_more = func.import_data(ids_more_data, cols_to_remove)

    if r == "yes":
        x_split1 = func.get_average_model(x_split1, ranges)
        x_split2 = func.get_average_model(x_split2, ranges)
        x_split3 = func.get_average_model(x_split3, ranges)
        x_split4 = func.get_average_model(x_split4, ranges)
        x_split5 = func.get_average_model(x_split5, ranges)
        x_more = func.get_average_model(x_more, ranges)

    if more == "yes":
        x_train = x_split1 + x_split2 + x_split3 + x_split4 + x_more
        y_train = y_split1 + y_split2 + y_split3 + y_split4 + y_more
        train_length = length_split1 + length_split2 + length_split3 + length_split4 + length_more
    else:
        x_train = x_split1 + x_split2 + x_split3 + x_split4
        y_train = y_split1 + y_split2 + y_split3 + y_split4
        train_length = length_split1 + length_split2 + length_split3 + length_split4

    x_cross = x_split5
    y_cross = y_split5
    cross_length = length_split5

    return x_train, y_train, train_length, x_cross, y_cross, cross_length


def train_model(x_train, y_train):
    sm = SMOTE(random_state=42)
    x_train, y_train = sm.fit_sample(x_train, y_train)

    classifier = MLPClassifier(hidden_layer_sizes=(200,), alpha=0.0001, random_state=1, tol=0.0000001)
    classifier.fit(x_train, y_train)

    return classifier


def model_accuracy(classifier, x_cross, y_cross, cross_length):
    prediction = classifier.predict(x_cross)

    # split prediction in per protein prediction
    i = 0
    performances = []
    for l in cross_length:
        pp = prediction[i:l]
        y_pp = y_cross[i:l]

        prec = metrics.precision_score(y_pp, pp, average=None)
        sensi = metrics.recall_score(y_pp, pp, average=None)
        f1 = metrics.f1_score(y_pp, pp, average=None)
        acc = metrics.accuracy_score(y_pp, pp)

        performance = [round(prec[0], 3), round(prec[1], 3), round(sensi[0], 3), round(sensi[1], 3), round(f1[0], 3),
                       round(f1[1], 3), round(acc, 3)]
        performances.append(performance)
        i = i + l

    mean_performances = np.mean(performances, axis=0)
    sd_performances = np.std(performances, axis=0)
    sd_performances = sd_performances / sqrt(999)

    mean_line = ""
    sd_line = ""
    for x in mean_performances:
        mean_line = mean_line + "\t" + str(x)
    mean_line = mean_line + "\t" + str(classifier.n_iter_) + "\t" + str(len(classifier.coefs_[0])) + "\n"

    for x in sd_performances:
        sd_line = sd_line + "\t" + str(x)
    sd_line = sd_line + "\n"

    return mean_line, sd_line


def prec_cov_curve(classifier, x_cross, y_cross, cross_length):
    proba = classifier.predict_proba(x_cross)

    cutoffs = range(0, 10000, 1)
    precision = []
    coverage = []

    for cut in cutoffs:
        # print(cut)
        float_cut = cut / 10000

        i = 0
        tmp_prec = []
        tmp_cov = []
        for l in cross_length:
            pp = proba[i:l]
            y_pp = y_cross[i:l]

            prediction = []
            for el in pp:
                if el[1] >= float_cut:
                    prediction.append(1)
                else:
                    prediction.append(0)
            prec = metrics.precision_score(y_pp, prediction, average=None)[1]
            cov = metrics.recall_score(y_pp, prediction, average=None)[1]
            tmp_prec.append(prec)
            tmp_cov.append(cov)

            i = i + 1

        mean_prec = np.mean(tmp_prec)
        mean_cov = np.mean(tmp_cov)

        precision.append(mean_prec)
        coverage.append(mean_cov)

    result = [coverage, precision]
    return result


def test_ovetraining(model, x_train, y_train, x_cross, y_cross):
    hidden_layers = ((10,), (50,), (100,), (200,), (300,), (500,), (800,), (1000,))
    # hidden_layers = ((1,), (700,))
    iterations = (20, 40, 60, 80, 100, 120, 140, 160, 300, 500)

    for layers in hidden_layers:
        print(layers)

        train_scores = []
        test_scores = []

        for iter in iterations:
            print(iter)
            classifier = MLPClassifier(hidden_layer_sizes=layers, max_iter=iter, tol=-100)
            classifier.fit(x_train, y_train)

            train_score = classifier.score(x_train, y_train)
            test_score = classifier.score(x_cross, y_cross)

            train_scores.append(train_score)
            test_scores.append(test_score)

            # print(train_score)
            # print(test_score)

            print(classifier.n_iter_)

        # plot figure

        fig = plt.figure()

        plt.title("Validation Curve for " + model + " " + str(layers))
        plt.xlabel("Iterations")
        plt.ylabel("Accuracy")
        plt.ylim(0.6, 1.0)
        lw = 2
        plt.plot(iterations, train_scores, label="Training score",
                 color="red", lw=lw)
        plt.plot(iterations, test_scores, label="Cross-validation score",
                 color="navy", lw=lw)
        plt.legend(loc="best")

        fig.savefig(
            "D:/Dropbox/masterthesis/thesis/plots/machine_learning/hidden_units/" + model + "_" + str(layers) + ".png")
        plt.close("all")
