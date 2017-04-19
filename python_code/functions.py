import csv
from sklearn import metrics, utils


def classifier_estimation(classifier, x_test, y_test):
    print()
    print('Test set has ', len([i for i in y_test if i is 1]), 'positives and ',
          len([i for i in y_test if i is 0]), ' negatives')
    print()
    print('Best classifier score')
    print()
    print(metrics.classification_report(y_test, classifier.predict(x_test), target_names=['non-binding', 'binding']))


def import_data(ids, cols_to_remove):
    positives = []
    negatives = []

    for id in ids:
        print(id)
        with open("data/single_files/" + id + ".pos") as p:
            for row in csv.reader(p, delimiter="\t"):
                # row = [i for j, i in enumerate(row) if j not in cols_to_remove]
                row = [float(i) for i in row]
                positives.append(row)
        with open("data/single_files/" + id + ".neg") as n:
            for row in csv.reader(n, delimiter="\t"):
                # row = [i for j, i in enumerate(row) if j not in cols_to_remove]
                row = [float(i) for i in row]
                negatives.append(row)

    data = positives + negatives
    labels = [1] * len(positives) + [0] * len(negatives)

    return data, labels


def calculate_performance(classifier, x_test, y_test):
    prec = metrics.precision_score(y_test, classifier.predict(x_test), average=None)
    sensi = metrics.recall_score(y_test, classifier.predict(x_test), average=None)
    f1 = metrics.f1_score(y_test, classifier.predict(x_test), average=None)
    acc = metrics.accuracy_score(y_test, classifier.predict(x_test))

    performance = [round(prec[0], 3), round(prec[1], 3), round(sensi[0], 3), round(sensi[1], 3), round(f1[0], 3),
                   round(f1[1], 3), round(acc, 3)]
    return performance


def bootstrapping(classifier, x_test, y_test):
    size = int(round(len(x_test) / 2))
    performances = []

    for i in range(0, 999):
        part_x, part_y = utils.resample(x_test, y_test, replace=True, n_samples=size)
        performance = calculate_performance(classifier, part_x, part_y)
        # string_performance = "\t".join([str(i) for i in tmp_performance])
        performances.append(performance)

    return performances
