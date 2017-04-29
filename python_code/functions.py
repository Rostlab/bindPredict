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
    data = []
    labels = []
    lengths = []

    for id in ids:
        length_pos = 0
        length_neg = 0

        with open("data/single_files/" + id + ".pos") as p:
            for row in csv.reader(p, delimiter="\t"):
                row = [i for j, i in enumerate(row) if j not in cols_to_remove]
                row = [float(i) for i in row]
                data.append(row)
                length_pos = length_pos + 1
        with open("data/single_files/" + id + ".neg") as n:
            for row in csv.reader(n, delimiter="\t"):
                row = [i for j, i in enumerate(row) if j not in cols_to_remove]
                row = [float(i) for i in row]
                data.append(row)
                length_neg = length_neg + 1

        tmp_labels = [1]*length_pos + [0]*length_neg
        labels.append(tmp_labels)
        length = length_pos+length_neg
        lengths.append(length)

    return data, labels, lengths


def get_average_model(data, ranges):
    new_data = []

    for i in range(0, len(data)):
        row = data[i]
        new_row = []
        for j in range(0, len(ranges)):
            curr_range = ranges[j]
            first = curr_range[0]
            last = curr_range[1] + 1
            avg = sum(row[first:last]) / len(row[first:last])
            new_row.append(avg)
        new_data.append(new_row)

    return new_data


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
