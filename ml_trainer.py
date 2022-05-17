from sklearn.model_selection import PredefinedSplit
from collections import defaultdict
import itertools as it

import torch

from architectures import CNN2Layers
from data_preparation import MyDataset
from assess_performance import ModelPerformance, PerformanceEpochs
from config import FileManager, GeneralInformation
from pytorchtools import EarlyStopping


class MLTrainer(object):

    def __init__(self, pos_weights, batch_size=406):
        self.batch_size = batch_size
        if torch.cuda.is_available():
            self.device = 'cuda:0'
        else:
            self.device = 'cpu'

        self.pos_weights = torch.tensor(pos_weights).to(self.device)

    def train_validate(self, params, train_ids, validation_ids, sequences, embeddings, labels, max_length,
                       verbose=True):
        """
        Train & validate predictor for one set of parameters and ids
        :param params:
        :param train_ids:
        :param validation_ids:
        :param sequences:
        :param embeddings:
        :param labels:
        :param max_length:
        :param verbose:
        :return:
        """

        model, train_performance, val_performance = self._train_validate(params, train_ids, validation_ids, sequences,
                                                                         embeddings, labels, max_length,
                                                                         verbose=verbose)

        train_loss, train_acc, train_prec, train_recall, train_f1, train_mcc = \
            train_performance.get_performance_last_epoch()
        val_loss, val_acc, val_prec, val_recall, val_f1, val_mcc = val_performance.get_performance_last_epoch()

        print("Train loss: {:.3f}, Prec: {:.3f}, Recall: {:.3f}, F1: {:.3f}, MCC: {:.3f}".format(train_loss,
                                                                                                 train_prec,
                                                                                                 train_recall,
                                                                                                 train_f1,
                                                                                                 train_mcc))
        print("Val loss: {:.3f}, Prec: {:.3f}, Recall: {:.3f}, F1: {:.3f}, MCC: {:.3f}".format(val_loss,
                                                                                               val_prec,
                                                                                               val_recall, val_f1,
                                                                                               val_mcc))

        return model

    def cross_validate(self, params, ids, fold_array, sequences, embeddings, labels, max_length, result_file):
        """
        Perform cross-validation to optimize hyperparameters
        :param params:
        :param ids:
        :param fold_array:
        :param sequences:
        :param embeddings:
        :param labels:
        :param max_length:
        :param result_file:
        :return:
        """
        ps = PredefinedSplit(fold_array)

        # create parameter grid
        param_sets = defaultdict(dict)
        sorted_keys = sorted(params.keys())
        param_combos = it.product(*(params[s] for s in sorted_keys))
        counter = 0
        for p in list(param_combos):
            curr_params = list(p)
            param_dict = dict(zip(sorted_keys, curr_params))
            param_sets[counter] = param_dict

            counter += 1

        best_score = 0
        best_params = dict()  # save best parameter set
        best_classifier = None  # save best classifier
        performance = defaultdict(list)  # save performance for each parameter combination

        params_counter = 1

        for p in param_sets.keys():
            curr_params = param_sets[p]
            print('{}\t{}'.format(params_counter, curr_params))

            model = None

            train_model_performance = ModelPerformance()
            val_model_performance = ModelPerformance()

            for train_index, test_index in ps.split():
                train_ids, validation_ids = ids[train_index], ids[test_index]

                model, train_performance, val_performance = self._train_validate(curr_params, train_ids, validation_ids,
                                                                                 sequences, embeddings, labels,
                                                                                 max_length)

                train_loss, train_acc, train_prec, train_recall, train_f1, train_mcc = \
                    train_performance.get_performance_last_epoch()
                val_loss, val_acc, val_prec, val_recall, val_f1, val_mcc = val_performance.get_performance_last_epoch()

                train_model_performance.add_single_performance(train_loss, train_acc, train_prec, train_recall,
                                                               train_f1, train_mcc)

                val_model_performance.add_single_performance(val_loss, val_acc, val_prec, val_recall, val_f1, val_mcc)

            # take average over all splits
            train_loss, train_acc, train_prec, train_recall, train_f1, train_mcc = \
                train_model_performance.get_mean_performance()
            val_loss, val_acc, val_prec, val_recall, val_f1, val_mcc = val_model_performance.get_mean_performance()

            performance['train_precision'].append(train_prec)
            performance['train_recall'].append(train_recall)
            performance['train_f1'].append(train_f1)
            performance['train_mcc'].append(train_mcc)
            performance['train_acc'].append(train_acc)
            performance['train_loss'].append(train_loss)

            performance['val_precision'].append(val_prec)
            performance['val_recall'].append(val_recall)
            performance['val_f1'].append(val_f1)
            performance['val_mcc'].append(val_mcc)
            performance['val_acc'].append(val_acc)
            performance['val_loss'].append(val_loss)

            for param in curr_params.keys():
                performance[param].append(curr_params[param])

            if val_f1 > best_score:
                best_score = val_f1
                best_params = curr_params
                best_classifier = model

            params_counter += 1

        FileManager.save_cv_results(performance, result_file)

        print(best_score)
        print(best_params)

        return best_classifier

    def _train_validate(self, params, train_ids, validation_ids, sequences, embeddings, labels, max_length,
                        verbose=True):
        """
        Train and validate bindEmbed21DL model
        :param params:
        :param train_ids:
        :param validation_ids:
        :param sequences:
        :param labels:
        :param max_length:
        :param verbose:
        :return:
        """

        # define data sets
        train_set = MyDataset(train_ids, embeddings, sequences, labels, max_length)
        validation_set = MyDataset(validation_ids, embeddings, sequences, labels, max_length)

        train_loader = torch.utils.data.DataLoader(train_set, batch_size=self.batch_size, shuffle=True, pin_memory=True,
                                                   worker_init_fn=GeneralInformation.seed_worker)
        validation_loader = torch.utils.data.DataLoader(validation_set, batch_size=self.batch_size, shuffle=True,
                                                        pin_memory=True, worker_init_fn=GeneralInformation.seed_worker)

        pos_weights = self.pos_weights.expand(max_length, 3)
        pos_weights = pos_weights.t()

        loss_fun = torch.nn.BCEWithLogitsLoss(reduction='none', pos_weight=pos_weights)
        sigm = torch.nn.Sigmoid()

        padding = int((params['kernel'] - 1) / 2)
        model = CNN2Layers(train_set.get_input_dimensions(), params['features'], params['kernel'], params['stride'],
                           padding, params['dropout'])
        model.to(self.device)

        optim_args = {'lr': params['lr'], 'betas': params['betas'], 'eps': params['eps'],
                      'weight_decay': params['weight_decay']}
        optimizer = torch.optim.Adamax(model.parameters(), **optim_args)

        checkpoint_file = 'checkpoint_early_stopping.pt'
        early_stopping = EarlyStopping(patience=10, delta=0.01, checkpoint_file=checkpoint_file)

        train_performance = PerformanceEpochs()
        validation_performance = PerformanceEpochs()

        num_epochs = 0

        for epoch in range(params['epochs']):
            if verbose:
                print("Epoch {}".format(epoch))

            train_loss = val_loss = 0
            train_loss_count = val_loss_count = 0
            train_tp = train_tn = train_fn = train_fp = 0
            val_tp = val_tn = val_fn = val_fp = 0

            train_acc = train_prec = train_rec = train_f1 = train_mcc = 0
            val_acc = val_prec = val_rec = val_f1 = val_mcc = 0

            # training
            model.train()
            for in_feature, target, loss_mask in train_loader:
                optimizer.zero_grad()
                in_feature = in_feature.to(self.device)
                in_feature_1024 = in_feature[:, :-1, :]

                target = target.to(self.device)
                loss_mask = loss_mask.to(self.device)
                pred = model.forward(in_feature_1024)

                # don't consider padded positions for loss calculation
                loss_el = loss_fun(pred, target)
                loss_el_masked = loss_el * loss_mask

                loss_norm = torch.sum(loss_el_masked)

                train_loss += loss_norm.item()
                train_loss_count += in_feature.shape[0]

                for idx, i in enumerate(in_feature):  # remove padded positions to calculate tp, fp, tn, fn
                    pred_i, target_i = GeneralInformation.remove_padded_positions(pred[idx], target[idx], i)

                    pred_i = sigm(pred_i)
                    tp, fp, tn, fn, acc, prec, rec, f1, mcc = \
                        train_performance.get_performance_batch(pred_i.detach().cpu(), target_i.detach().cpu())
                    train_tp += tp
                    train_fp += fp
                    train_tn += tn
                    train_fn += fn

                    train_acc += acc
                    train_prec += prec
                    train_rec += rec
                    train_f1 += f1
                    train_mcc += mcc

                loss_norm.backward()
                optimizer.step()

            # validation
            model.eval()
            with torch.no_grad():
                for in_feature, target, loss_mask in validation_loader:
                    in_feature = in_feature.to(self.device)
                    in_feature_1024 = in_feature[:, :-1, :]

                    target = target.to(self.device)
                    loss_mask = loss_mask.to(self.device)

                    pred = model.forward(in_feature_1024)

                    # don't consider padded position for loss calculation
                    loss_el = loss_fun(pred, target)
                    loss_el_masked = loss_el * loss_mask
                    val_loss += torch.sum(loss_el_masked).item()
                    val_loss_count += in_feature.shape[0]

                    for idx, i in enumerate(in_feature):  # remove padded positions to calculate tp, fp, tn, fn
                        pred_i, target_i = GeneralInformation.remove_padded_positions(pred[idx], target[idx], i)

                        pred_i = sigm(pred_i)

                        tp, fp, tn, fn, acc, prec, rec, f1, mcc = \
                            train_performance.get_performance_batch(pred_i.detach().cpu(), target_i.detach().cpu())
                        val_tp += tp
                        val_fp += fp
                        val_tn += tn
                        val_fn += fn

                        val_acc += acc
                        val_prec += prec
                        val_rec += rec
                        val_f1 += f1
                        val_mcc += mcc

            train_loss = train_loss / (train_loss_count * 3)
            val_loss = val_loss / (val_loss_count * 3)

            train_acc = train_acc / train_loss_count
            train_prec = train_prec / train_loss_count
            train_rec = train_rec / train_loss_count
            train_f1 = train_f1 / train_loss_count
            train_mcc = train_mcc / train_loss_count

            val_acc = val_acc / val_loss_count
            val_prec = val_prec / val_loss_count
            val_rec = val_rec / val_loss_count
            val_f1 = val_f1 / val_loss_count
            val_mcc = val_mcc / val_loss_count

            if verbose:
                print("Train loss: {:.3f}, Prec: {:.3f}, Recall: {:.3f}, F1: {:.3f}, MCC: {:.3f}".format(train_loss,
                                                                                                         train_prec,
                                                                                                         train_rec,
                                                                                                         train_f1,
                                                                                                         train_mcc))
                print('TP: {}, FP: {}, TN: {}, FN: {}'.format(train_tp, train_fp, train_tn, train_fn))
                print("Val loss: {:.3f}, Prec: {:.3f}, Recall: {:.3f}, F1: {:.3f}, MCC: {:.3f}".format(val_loss,
                                                                                                       val_prec,
                                                                                                       val_rec, val_f1,
                                                                                                       val_mcc))
                print('TP: {}, FP: {}, TN: {}, FN: {}'.format(val_tp, val_fp, val_tn, val_fn))

            # append average performance for this epoch
            train_performance.add_performance_epoch(train_loss, train_mcc, train_prec, train_rec, train_f1, train_acc)
            validation_performance.add_performance_epoch(val_loss, val_mcc, val_prec, val_rec, val_f1, val_acc)

            num_epochs += 1

            # stop training if F1 score doesn't improve anymore
            if 'early_stopping' in params.keys() and params['early_stopping']:
                eval_val = val_f1 * (-1)
                # eval_val = val_loss
                early_stopping(eval_val, model, verbose)
                if early_stopping.early_stop:
                    break

        model = torch.load(checkpoint_file)

        return model, train_performance, validation_performance
