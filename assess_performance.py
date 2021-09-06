import torch
import math
import numpy as np
from scipy.stats import t


class ModelPerformance(object):

    """Wrapper to save performances values per model"""

    def __init__(self):

        self.losses = []

        self.precisions = []
        self.recalls = []
        self.f1s = []

        self.neg_precisions = []
        self.neg_recalls = []
        self.neg_f1s = []

        self.accs = []
        self.mccs = []

        self.confusion_matrix = {'tp': 0, 'fp': 0, 'tn': 0, 'fn': 0}
        self.cross_prediction = [0, 0, 0, 0]

        self.coverage = 0
        self.predictions = {'bound_predicted': 0, 'bound_not_predicted': 0, 'not_bound_predicted': 0,
                            'not_bound_not_predicted': 0}

    def add_single_performance(self, loss, acc, prec, recall, f1, mcc):
        """Add performance for one protein"""

        self.losses.append(loss)
        self.accs.append(acc)
        self.precisions.append(prec)
        self.recalls.append(recall)
        self.f1s.append(f1)
        self.mccs.append(mcc)

    def add_single_performance_negatives(self, loss, acc, prec, recall, f1, neg_prec, neg_recall, neg_f1, mcc):
        """Add performance for one protein (incl. negative performance values)"""

        self.add_single_performance(loss, acc, prec, recall, f1, mcc)
        self.neg_precisions.append(neg_prec)
        self.neg_recalls.append(neg_recall)
        self.neg_f1s.append(neg_f1)

    def add_confusion_matrix(self, tp, fp, tn, fn):
        """Add TP, FP, TN,FN to confusion matrix"""

        self.confusion_matrix['tp'] += tp
        self.confusion_matrix['fp'] += fp
        self.confusion_matrix['tn'] += tn
        self.confusion_matrix['fn'] += fn

        if (tp + fp) > 0:
            self.coverage += 1
            if (tp + fn) > 0:
                self.predictions['bound_predicted'] += 1
            else:
                self.predictions['not_bound_predicted'] += 1
        else:
            if fn > 0:
                self.predictions['bound_not_predicted'] += 1
            else:
                self.predictions['not_bound_not_predicted'] += 1

    def add_cross_prediction(self, cross_prediction):
        """Add cross prediction"""
        for i in range(0, len(cross_prediction)):
            self.cross_prediction[i] += cross_prediction[i]

    def get_mean_performance(self):
        """Calculate average performance values"""
        loss = np.average(self.losses)
        acc = np.average(self.accs)
        precision = np.average(self.precisions)
        recall = np.average(self.recalls)
        f1 = np.average(self.f1s)
        mcc = np.average(self.mccs)

        return loss, acc, precision, recall, f1, mcc

    def get_mean_ci_performance(self):
        """Calculate average performance values and 95% CIs"""
        acc, acc_ci = ModelPerformance._get_mean_ci(self.accs)
        recall, recall_ci = ModelPerformance._get_mean_ci(self.recalls)
        prec, prec_ci = ModelPerformance._get_mean_ci(self.precisions)
        f1, f1_ci = ModelPerformance._get_mean_ci(self.f1s)
        mcc, mcc_ci = ModelPerformance._get_mean_ci(self.mccs)

        return acc, prec, recall, f1, mcc, acc_ci, prec_ci, recall_ci, f1_ci, mcc_ci

    def get_mean_ci_performance_negatives(self):
        """Calculate average performance values for negative class and 95% CIs"""

        neg_recall, neg_recall_ci = ModelPerformance._get_mean_ci(self.neg_recalls)
        neg_prec, neg_prec_ci = ModelPerformance._get_mean_ci(self.neg_precisions)
        neg_f1, neg_f1_ci = ModelPerformance._get_mean_ci(self.neg_f1s)

        return neg_prec, neg_recall, neg_f1, neg_prec_ci, neg_recall_ci, neg_f1_ci

    @staticmethod
    def _get_mean_ci(vec):
        """
        Calculate mean and 95% CI for a given vector
        :param vec: vector
        :return: mean and ci
        """
        mean = round(np.average(vec), 3)
        if len(vec) > 1:
            ci = round(np.std(vec)/math.sqrt(len(vec)) * t.ppf((1 + 0.95) / 2, len(vec)), 3)
        else:
            ci = 0

        return mean, ci


class PerformanceEpochs(object):
    """
    Wrapper to save performance values per epoch
    """

    def __init__(self):

        self.loss_epochs = []
        self.mcc_epochs = []
        self.prec_epochs = []
        self.recall_epochs = []
        self.f1_epochs = []
        self.acc_epochs = []

    def get_performance_last_epoch(self):
        """Get performance for last epoch"""
        loss = self.loss_epochs[-1]
        mcc = self.mcc_epochs[-1]
        prec = self.prec_epochs[-1]
        recall = self.recall_epochs[-1]
        f1 = self.f1_epochs[-1]
        acc = self.acc_epochs[-1]

        return loss, acc, prec, recall, f1, mcc

    def add_performance_epoch(self, loss, mcc, prec, recall, f1, acc):
        """Add performance for one epoch"""
        self.loss_epochs.append(loss)
        self.acc_epochs.append(acc)
        self.mcc_epochs.append(mcc)
        self.f1_epochs.append(f1)
        self.prec_epochs.append(prec)
        self.recall_epochs.append(recall)

    @staticmethod
    def get_performance_batch(pred, target):
        """Calculate performance for one batch"""

        tp, fp, tn, fn = PerformanceAssessment.evaluate_per_residue_torch(pred, target)
        acc, prec, rec, f1, mcc = PerformanceAssessment.calc_performance_measurements(tp, fp, tn, fn)

        return tp, fp, tn, fn, acc, prec, rec, f1, mcc


class PerformanceAssessment(object):

    @staticmethod
    def evaluate_per_residue_torch(prediction, target):
        """Calculate tp, fp, tn, fn for tensor"""
        # reduce prediction & target to one dimension
        prediction = prediction.t()
        target = target.t()
        prediction = torch.sum(torch.ge(prediction, 0.5), 1)
        target = torch.sum(torch.ge(target, 0.5), 1)

        # get confusion matrix
        tp = torch.sum(torch.ge(prediction, 0.5) * torch.ge(target, 0.5))
        tn = torch.sum(torch.lt(prediction, 0.5) * torch.lt(target, 0.5))
        fp = torch.sum(torch.ge(prediction, 0.5) * torch.lt(target, 0.5))
        fn = torch.sum(torch.lt(prediction, 0.5) * torch.ge(target, 0.5))

        return tp, fp, tn, fn

    @staticmethod
    def calc_performance_measurements(tp, fp, tn, fn):
        """Calculate precision, recall, f1, mcc, and accuracy"""

        tp = float(tp)
        fp = float(fp)
        fn = float(fn)
        tn = float(tn)

        recall = prec = f1 = mcc = 0
        acc = round((tp + tn) / (tp + tn + fn + fp), 3)

        if tp > 0 or fn > 0:
            recall = round(tp / (tp + fn), 3)
        if tp > 0 or fp > 0:
            prec = round(tp / (tp + fp), 3)
        if recall > 0 or prec > 0:
            f1 = round(2 * recall * prec / (recall + prec), 3)
        if (tp > 0 or fp > 0) and (tp > 0 or fn > 0) and (tn > 0 or fp > 0) and (tn > 0 or fn > 0):
            mcc = round((tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)), 3)

        return acc, prec, recall, f1, mcc

    @staticmethod
    def combine_protein_performance(proteins, cutoff, labels):
        """
        Calculate final model performance from per-protein performances
        :param proteins:
        :param cutoff:
        :param labels:
        :return:
        """
        model_performances = {'overall': ModelPerformance(), 'metal': ModelPerformance(),
                              'nucleic': ModelPerformance(), 'small': ModelPerformance()}

        for p in proteins.keys():
            prot = proteins[p]
            prot.set_labels(labels[p])

            # calculate performance for protein
            performance = prot.calc_performance_measurements(cutoff)

            for k in performance.keys():
                model_performance = model_performances[k]
                if 'tp' in performance[k].keys():
                    tp = performance[k]['tp']
                    fp = performance[k]['fp']
                    fn = performance[k]['fn']
                    tn = performance[k]['tn']
                    if k == 'overall' and (tp + fn) == 0:
                        print('No residues annotated as binding for {}'.format(p))
                    if (tp + fp + fn) > 0:
                        # only add performance if this protein binds or if one residue was predicted to bind
                        model_performance.add_single_performance(0, performance[k]['acc'], performance[k]['prec'],
                                                                 performance[k]['recall'], performance[k]['f1'],
                                                                 performance[k]['mcc'])
                        model_performance.add_confusion_matrix(tp, fp, tn, fn)
                    else:
                        model_performance.predictions['not_bound_not_predicted'] += 1

        return model_performances

    @staticmethod
    def write_performance_results(model_performances, out_file):
        """Write average performance"""
        with open(out_file, 'w') as out:
            out.write('Type\ttp\tfp\ttn\tfn\tprec\tprec.ci\trecall\trecall.ci\tf1\tf1.ci\tmcc\tmcc.ci\tacc\tacc.ci\n')
            for k in model_performances.keys():
                model_performance = model_performances[k]
                acc, pr, rec, f1, mcc, acc_ci, pr_ci, rec_ci, f1_ci, mcc_ci = \
                    model_performance.get_mean_ci_performance()

                confusion_matrix = model_performance.confusion_matrix
                tp = confusion_matrix['tp']
                fp = confusion_matrix['fp']
                tn = confusion_matrix['tn']
                fn = confusion_matrix['fn']

                out.write('{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}'
                          '\t{:.3f}\n'.format(k, tp, fp, tn, fn, pr, pr_ci, rec, rec_ci, f1, f1_ci, mcc, mcc_ci, acc,
                                              acc_ci))

    @staticmethod
    def print_performance_results(model_performances):
        """Print average performance"""

        for k in model_performances.keys():
            print(k)
            model_performance = model_performances[k]
            acc, pr, rec, f1, mcc, acc_ci, pr_ci, rec_ci, f1_ci, mcc_ci = \
                model_performance.get_mean_ci_performance()

            cov_proteins = model_performance.coverage
            if len(model_performance.accs) > 0:
                cov_percentage = cov_proteins / len(model_performance.accs)
            else:
                cov_percentage = 0.0

            confusion_matrix = model_performance.confusion_matrix
            predictions = model_performance.predictions

            print('CovOneBind: {} ({:.3f})'.format(cov_proteins, cov_percentage))
            print('Bound: With predictions: {}, Without predictions: {}\nNot Bound: With predictions: {}, '
                  'Without predictions: {}'.format(predictions['bound_predicted'], predictions['bound_not_predicted'],
                                                   predictions['not_bound_predicted'],
                                                   predictions['not_bound_not_predicted']))
            print('TP: {}, FP: {}, TN: {}, FN: {}'.format(confusion_matrix['tp'], confusion_matrix['fp'],
                                                          confusion_matrix['tn'], confusion_matrix['fn']))
            print("Prec: {:.3f} +/- {:.3f}, Recall: {:.3f} +/- {:.3f}, F1: {:.3f} +/- {:.3f}, "
                  "MCC: {:.3f} +/- {:.3f}, Acc: {:.3f} +/- {:.3f}".format(pr, pr_ci, rec, rec_ci, f1, f1_ci, mcc,
                                                                          mcc_ci, acc, acc_ci))

    @staticmethod
    def print_cross_prediction_results(model_performances):
        """Print cross-predictions"""

        for k in model_performances.keys():
            print(k)
            cross_predictions = model_performances[k].cross_prediction
            print('Metal: {}, Nucleic: {}, Small: {}, Non-Binding: {}'.format(cross_predictions[0],
                                                                              cross_predictions[1],
                                                                              cross_predictions[2],
                                                                              cross_predictions[3]))
