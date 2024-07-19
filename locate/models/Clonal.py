import itertools

import torch
import torch.nn as nn
from torch.distributions import constraints

from locate.models.Model import Model
from locate.likelihoods import ClonalLikelihood

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.infer import infer_discrete, config_enumerate
from pyro.infer.autoguide import AutoDelta
from pyro.ops.indexing import Vindex
from pyro.util import ignore_jit_warnings


# Here I created a new class that inherits from dict, it ha an unsqueeze method that unsqueezes all the tensors in the dict 
# It is actually more of a dictionary generator than a proper dictionary

class SqueezableDict(dict):
    
    def __init__(self, *args, **kwargs):
        super(SqueezableDict, self).__init__(*args, **kwargs)
        self.data = [v for v in self.values() if v is not None ][0]
        self.dim = [v for v in self.values() if v is not None ][0].dim
        self.shape = [v for v in self.values() if v is not None ][0].shape
    
    def unsqueeze(self, dim):
        return {k:v.unsqueeze(dim) for k,v in self.items() if v is not None }
    


class Clonal(Model):

    # parameters of the HMM + other stuff specific to CNA calling

    params = {'jumping_prob' : 1e-2,                                              
              'init_probs': torch.tensor([0.2, 0.3, 0.3, 0.1, 0.1]), 
              'hidden_dim': 5, 
              'atac_lk' : "P", 
              "CUDA" : False, 
              "multiome" : False, 
              "prior_ploidy" : 2, 
              "prior_purity" : 1, 
              "fix_BAF_emissions" : False, 
              "fix_purity" : True, 
              "N_bins" : 3000, 
              "scaling_factors" : torch.tensor([1.,1.,1.,1.])}
    
    # This has to  be generalized and left empty 
    data_name = set(['exp_rna', 'cov_atac', 'baf_rna', 'baf_atac'])

    
    def __init__(self, data_dict):
        self._params = self.params.copy()
        self._data = None
        self._name = "SegmentAnything"
        super().__init__(data_dict, self.data_name)
        # internal structure is a SqueezableDict
        self._data = SqueezableDict({k:v for k,v in self._data.items()})
            
    def model(self, i = 1,  *args, **kwargs):
        
        n_sequences, length = 0, 0
        tot = 0

        # check what data types we have, all the emissions are hardcoded for now
        if ((self._data["exp_rna"]  is not None  or self._data["cov_atac"]  is not None )) \
        and ((self._data["baf_rna"]   is not None or self._data["baf_atac"]  is not None)):
            Major, minor, tot, x = self.get_Major_minor()

            # if we have the baf go for an allele specific model, otherwise encode simply the number of chromosomes
            self._params["has_baf"] = True
            self._params["Major"] = Major
            self._params["minor"] = minor
            probs_x = pyro.sample(
                "probs_x",
                dist.Dirichlet((1 - self._params["jumping_prob"]) * torch.eye(x.shape[0]) + self._params["jumping_prob"]).to_event(1),
            )
            
        else:
            self._params["has_baf"] = False
            probs_x = pyro.sample(
                "probs_x",
                dist.Dirichlet((1 - self._params["jumping_prob"]) * torch.eye(self._params["hidden_dim"]) + self._params["jumping_prob"]).to_event(1),
            )
            x = torch.arange(0, self._params["hidden_dim"]).long()
            tot = x + 1
            minor, Major = None, None

        # From here on we have the same code for the two models
        # Hardcoded emission parameters

        if self._data["exp_rna"] is not None:
            length, n_sequences = self._data["exp_rna"].shape
            exp_rna_mu = pyro.sample(
                "rna_mu",
                dist.Normal(0., 1.).expand([x.shape[0]]).to_event(1),
            )
            exp_rna_sd = pyro.sample(
                "rna_sd",
                dist.InverseGamma(1., 1.)
            )
        else:
            exp_rna_mu, exp_rna_sd = None, None
        
        if self._data["cov_atac"] is not None:
            length, n_sequences  = self._data["cov_atac"].shape
            cov_atac_mu = pyro.sample(
                "atac_mu",
                dist.Normal(torch.mean(self._data["cov_atac"]), torch.std(self._data["cov_atac"])),
            )
            if self._params["atac_lk"] == "NB":
                cov_atac_overdispersion = pyro.sample(
                "atac_overdispersion",
                dist.Uniform(torch.mean(self._data["cov_atac"]), torch.sd(self._data["cov_atac"])).expand([x.shape[0]]).to_event(1),
            )
            else:
                cov_atac_overdispersion = None
        else:
            cov_atac_mu = None
            
        if self._data["baf_rna"] is not None:
            length, n_sequences  = self._data["baf_rna"].shape
            baf_rna_number_of_trials = pyro.sample(
                "baf_rna_number_of_trials",
                dist.Uniform(1, 10000),
            ) 
        else:
            baf_rna_number_of_trials = None
            
        if self._data["baf_atac"] is not None:
            length, n_sequences  = self._data["baf_atac"].shape
            baf_atac_number_of_trials = pyro.sample(
                "baf_atac_number_of_trials",
                dist.Uniform(1, 10000),
            )
        else:
            baf_atac_number_of_trials = None
            
        if self._params["fix_purity"]:
            purity = self._params["prior_purity"]
        else:
            purity = pyro.sample("purity", dist.Uniform(0.,1.))
            
        self._params["N_bins"] = length
                
        with pyro.plate("sequences", n_sequences):
            
            init_logits = self._params["init_probs"].log()
            trans_logits = probs_x.log()
            
            # the actual likelihood
            ### WARNING!!! I think there is something wrong here!!! ###
            with ignore_jit_warnings():
                obs_dist = MultiomeLikelihood(
                 x = x.unsqueeze(-1),
                 cov_atac_mu = cov_atac_mu, 
                 cov_atac_overdispersion = cov_atac_overdispersion, 
                 exp_rna_mu = exp_rna_mu,
                 exp_rna_sd = exp_rna_sd,
                 Major = Major,
                 minor = minor,
                 tot = tot,
                 baf_rna_number_of_trials = baf_rna_number_of_trials,
                 baf_atac_number_of_trials = baf_atac_number_of_trials,
                 scaling_factors = self._params["scaling_factors"],
                 purity = purity, 
                 atak_lk = "P",
                 batch_shape = [x.shape[0],length]
                ).to_event(1)

                hmm_dist = dist.DiscreteHMM(init_logits, trans_logits, obs_dist)
            pyro.sample("y", hmm_dist, obs=self._data)
                
    # Autoguide, easy peasy lemon squeezy
    def guide(self, *args, **kwargs):
        return AutoDelta(poutine.block(self.model, hide_fn=lambda msg: msg["name"].startswith("x")))
    
    # If allele specific model map the states to the major and minor alleles
    def get_Major_minor(self):
        combinations = list(itertools.product(range(1, self._params["hidden_dim"]+1), range(self._params["hidden_dim"]+1)))
        # Separate the major and minor alleles
        major_alleles = [combination[0] for combination in combinations]
        minor_alleles = [combination[1] for combination in combinations]

        # Convert the lists into PyTorch tensors
        major_allele_tensor = torch.tensor(major_alleles).long()
        minor_allele_tensor = torch.tensor(minor_alleles).long()
        x = torch.tensor(list(range(len(major_alleles)))).long()
        tot = major_allele_tensor + minor_allele_tensor
        
        return major_allele_tensor, minor_allele_tensor, tot,  x
    
    
    # Model 2 is the same as model 1 but with enumeration, is it used to get MAP estimates of the states, becouse Pyro is a MESSED UP library
    @infer_discrete(first_available_dim=-2, temperature=0)
    @config_enumerate
    def model_2(self,learned_params):
        n_sequences, length = 0, 0
        minor, Major = None, None
        x = torch.arange(1., self._params["hidden_dim"] + 1)
        tot = x
        if ((self._data["exp_rna"]  is not None  or self._data["cov_atac"]  is not None )) \
        and ((self._data["baf_rna"]   is not None or self._data["baf_atac"]  is not None)):
            Major, minor, tot, x = self.get_Major_minor()
            
        probs_x = learned_params['probs_x']
            

            

        if self._data["exp_rna"] is not None:
            length, n_sequences = self._data["exp_rna"].shape
            exp_rna_mu = learned_params['rna_mu']
            exp_rna_sd = learned_params['rna_sd']
        else:
            exp_rna_mu, exp_rna_sd = None, None
        
        if self._data["cov_atac"] is not None:
            length, n_sequences  = self._data["cov_atac"].shape
            cov_atac_mu  = learned_params['atac_mu']
            if self._params["atac_lk"] == "NB":
                cov_atac_overdispersion = learned_params['atac_overdispersion']
            else:
                cov_atac_overdispersion = None
        else:
            cov_atac_mu = None
            
        if self._data["baf_rna"] is not None:
            length, n_sequences  = self._data["baf_rna"].shape
            baf_rna_number_of_trials = learned_params['baf_rna_number_of_trials']
        else:
            baf_rna_number_of_trials = None
            
        if self._data["baf_atac"] is not None:
            length, n_sequences  = self._data["baf_atac"].shape
            baf_atac_number_of_trials = learned_params['baf_atac_number_of_trials']
        else:
            baf_atac_number_of_trials = None
            
        if self._params["fix_purity"]:
            purity = self._params["prior_purity"]
        else:
            purity =  learned_params['purity']
            
                
        with pyro.plate("sequences", n_sequences, dim = -1):
            
            
            x = [0]
            for t in pyro.markov(range(length)):
                x_new = pyro.sample(
                "x_{}".format(t),
                dist.Categorical(probs_x[x[t]])
                )
                x.append(x_new)
                #print(x_new.shape)
                
                input_tmp = SqueezableDict({k:v[t,:] for k,v in self._data.items() if v is not None})
                
        
        
                if self._data["cov_atac"] is not None:
                    mean = purity * cov_atac_mu * (tot[x_new]) + (1 - purity) *cov_atac_mu * 2
                    if self._params["atac_lk"] == "NB":
                        cov_atac_lk = pyro.factor("y_atac_cov_{}".format(t), dist.NegativeBinomial(mean, cov_atac_overdispersion[x.squeeze(-1)]).log_prob(self._data["cov_atac"][t,:]))
                    else:
                        cov_atac_lk = pyro.factor("y_atac_cov_{}".format(t), dist.Poisson(mean).log_prob(self._data["cov_atac"][t,:]))
                if self._data["exp_rna"] is not None:
                    mean = purity * exp_rna_mu[x_new] + (1 - purity) * exp_rna_mu[1]
                    exp_rna_lk = pyro.factor("y_rna_exp_{}".format(t) , dist.Normal(mean, exp_rna_sd).log_prob(self._data["exp_rna"][t,:]))
                        

                if  self._data["baf_rna"] is not None:

                    prob_tum =  (minor[x_new]  / (Major[x_new] + minor[x_new]))  + 1e-6
                    prob = purity * prob_tum + 0.5 * (1 - purity)
                    baf_rna_lk = pyro.factor("y_rna_baf_{}".format(t) , dist.Beta(prob * baf_rna_number_of_trials, 
                                           (1 - prob) * baf_rna_number_of_trials).log_prob(self._data["baf_rna"][t,:]))



                if  self._data["baf_atac"] is not None:

                    prob_tum =  (minor[x_new]  / (Major[x_new] + minor[x_new])) + 1e-6
                    prob = purity * prob_tum + 0.5 * (1 - purity)
                    baf_atac_lk = pyro.factor("y_atac_baf_{}".format(t) , dist.Beta(prob * baf_atac_number_of_trials, 
                                            (1 - prob) * baf_atac_number_of_trials).log_prob(self._data["baf_atac"][t,:]))
            return x

           
        