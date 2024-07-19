
import pyro
import torch

from pyro.infer import SVI, infer_discrete
from pyro import poutine

from tqdm import trange

# from segmenter.building_blocks import *
# from segmenter.model_selection import *
# from locate.utils import retrieve_params
# from segmenter.stopping_criteria import all_stopping_criteria


class LOCATE:
    """ 
    The interface class
    
    Basically it takes a model, an optimizer and a loss function and provides a functions to run the inference and get the parameters
    back masking the differences between models.
    """
    
    def __init__(self,model = None, optimizer = None, loss = None, inf_type = SVI, CUDA = False):
        self._model_fun = model
        self._optimizer = optimizer
        self._loss = loss
        self._model = None
        self._inf_type = inf_type
        self._model_trained = None
        self._guide_trained = None
        self._loss_trained = None
        self._model_string = None
        self._CUDA = CUDA
        
        if self._CUDA:
            torch.set_default_tensor_type('torch.cuda.FloatTensor')
        else:
            torch.set_default_tensor_type('torch.FloatTensor')
        
        self.params_history = {}


    def __repr__(self):

        if self._model is None:
            dictionary = {"Model": self._model_fun,
                    "Optimizer": self._optimizer,
                    "Loss": self._loss,
                    "Inference strategy": self._inf_type
                    }
        else:
            dictionary = {"Model" : self._model_fun,
                    "Data" : self._model._data,
                    "Model params": self._model._params,
                    "Optimizer" :self._optimizer,
                    "Loss" : self._loss,
                    "Inference strategy" :  self._inf_type
                    }

        return "\n".join("{}:\t{}".format(k, v) for k, v in dictionary.items())

    def initialize_model(self, data):
        assert self._model_fun is not None
        self._model = self._model_fun(data)
        self._model_string = type(self._model).__name__

    def set_optimizer(self, optimizer):
        self._optimizer = optimizer

    def set_model(self, model):
        self._model_fun = model

    def set_loss(self, loss):
        self._loss = loss

    def set_model_params(self, param_dict):
        if self._CUDA:
            param_dict['CUDA'] = True
        self._model.set_params(param_dict)

    def run(self, steps,param_optimizer = {'lr' : 0.05}, e = 0.01, patience = 5, param_loss = None, seed = 3):

        """ This function runs the inference of non-categorical parameters

          This function performs a complete inference cycle for the given tuple(model, optimizer, loss, inference modality).
          For more info about the parameters for the loss and the optimizer look at `Optimization <http://docs.pyro.ai/en/stable/optimization.html>`_.
          and `Loss <http://docs.pyro.ai/en/stable/inference_algos.html>`_.

          Not all the the combinations Optimize-parameters and Loss-parameters have been tested, so something may
          not work (please open an issue on the GitHub page).


          Args:
              steps (int): Number of steps
              param_optimizer (dict):  A dictionary of paramaters:value for the optimizer
              param_loss (dict): A dictionary of paramaters:value for the loss function
              seed (int): seed to be passed to  pyro.set_rng_seed
              MAP (bool): Perform learn a Delta distribution over the outer layer of the model graph
              verbose(bool): show loss for each step, if false the functions just prints a progress bar
              BAF(torch.tensor): if provided use BAF penalization in the loss

          Returns:
              list: loss (divided by sample size) for each step


          """

        pyro.set_rng_seed(seed)
        pyro.clear_param_store()

        expose_vec = export_switch(self._model)
        model = self._model.model
        guide = self._model.guide(expose_vec)


        optim = self._optimizer(param_optimizer)
        elbo = self._loss(**param_loss) if param_loss is not None else self._loss()
        svi = self._inf_type(model, guide, optim, elbo)

        num_observations = 0
        num_segments = 0

        #if "data_rna" in self._model._data:
        #    num_observations += self._model._data['data_rna'].shape[1]
        #    num_segments += self._model._data['data_rna'].shape[0]

        #if "data_atac" in self._model._data:
        #    num_observations += self._model._data['data_atac'].shape[1]
        #    num_segments += self._model._data['data_atac'].shape[0]

        loss = [None] * steps

        #loss[0] = svi.step(i = 1) / (num_observations * num_segments)
        loss[0] = svi.step(i = 1)
        elb = loss[0]
        new_w = retrieve_params()

        t = trange(steps, desc='Bar desc', leave=True)

        conv = 0
        for step in t:
            t.set_description('ELBO: {:.9f}  '.format(elb))
            t.refresh()

            #loss[step] = svi.step(i = step + 1) / (num_observations * num_segments)
            loss[step] = svi.step(i = step + 1)
            elb = loss[step]

            old_w, new_w = new_w, retrieve_params()

            stop, diff_params = all_stopping_criteria(old_w, new_w, e, step)
            self.params_history.update({step : diff_params})

            if stop:
                if (conv == patience):
                    print('Reached convergence', flush = True)
                    break
                else:
                    conv = conv + 1
            else:
                conv = 0
        print("", flush = True)
        self._model_trained = model
        self._guide_trained = guide
        self._loss_trained = loss
        print("Done!", flush=True)
        return loss, num_observations


    def learned_parameters(self):

        """ Return all the estimated  parameter values

            Calls the right set of function for retrieving learned parameters according to the model type
            If posterior=True all the other parameters are just passed to :func:`~congas.core.Interface.inference_categorical_posterior`

            Args:

              posterior: learn posterior assignement (if false estimate MAP via Viterbi-like MAP inference)


            Returns:
              dict: parameter:value dictionary
        """


        params = self._guide_trained()
        
        if self._model._name != "SegmentAnything":

            if "data_rna" in self._model._data:
                if self._model._params["likelihood_rna"]  in ["G", "N"]:
                    params["segment_factor_rna"] = torch.ones(self._model._data['data_rna'].shape[0])

            if "data_atac" in self._model._data:
                if self._model._params["likelihood_atac"]  in ["G", "N"]:
                    params["segment_factor_atac"] = torch.ones(self._model._data['data_atac'].shape[0])

            params["CNA"] = torch.argmax(params["CNV_probabilities"] , axis=(2)) + 1

            print("", flush=True)

            print("Computing assignment probabilities", flush=True)
            discrete_params = self._model.calculate_cluster_assignements(params)
        else:
            print(params)
            map_states =  self._model.model_2(self._model, learned_params = params)
            if self._model._params["has_baf"]:
                discrete_params = {"CN_Major" : self._model._params["Major"][torch.tensor(map_states)[1:].long()],
                                  "CN_minor" : self._model._params["minor"][torch.tensor(map_states)[1:].long()]}
            else:
                discrete_params = {"CN_tot" : torch.tensor(map_states)[1:]}
            
        if self._CUDA:
            trained_params_dict = {i : params[i].cpu().detach().numpy() for i in params}
            discrete_params = {i : discrete_params[i].cpu().detach().numpy() for i in discrete_params}
        else:
            trained_params_dict = {i : params[i].detach().numpy() for i in params}
            discrete_params = {i : discrete_params[i].detach().numpy() for i in discrete_params}

        all_params =  {**trained_params_dict,**discrete_params}

        return all_params







