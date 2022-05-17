from data_preparation import MyDataset, ProteinResults
from config import GeneralInformation

import torch


class MLPredictor(object):

    def __init__(self, model):
        if torch.cuda.is_available():
            self.device = 'cuda:0'
        else:
            self.device = 'cpu'

        self.model = model.to(self.device)

    def predict_per_protein(self, ids, sequences, embeddings, labels, max_length):
        """
        Generate predictions for each protein from a list
        :param ids:
        :param sequences:
        :param embeddings:
        :param labels:
        :param max_length:
        :return:
        """

        validation_set = MyDataset(ids, sequences, embeddings, labels, max_length, protein_prediction=True)
        validation_loader = torch.utils.data.DataLoader(validation_set, batch_size=1, shuffle=True, pin_memory=True)
        sigm = torch.nn.Sigmoid()

        proteins = dict()
        for features, target, mask, prot_id in validation_loader:
            prot_id = prot_id[0]

            self.model.eval()
            with torch.no_grad():
                features = features.to(self.device)
                target = target.to(self.device)

                features_1024 = features[..., :-1, :]
                pred = self.model.forward(features_1024)
                pred = sigm(pred)

                pred = pred.squeeze()
                target = target.squeeze()
                features = features.squeeze()

                pred_i, target_i = GeneralInformation.remove_padded_positions(pred, target, features)
                pred_i = pred_i.detach().cpu()

                prot = ProteinResults(prot_id)
                prot.set_predictions(pred_i)
                proteins[prot_id] = prot

        return proteins
