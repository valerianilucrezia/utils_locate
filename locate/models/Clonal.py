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
              'init_probs': torch.tensor([0.1, 0.1, 0.1, 0.1, 0.1]), 
              'hidden_dim': 3, 
              "CUDA" : False, 
              "prior_ploidy" : 2, 
              "prior_purity" : 0.9, 
              "fix_purity" : True, 
              "fix_ploidy" : True, 
              "scaling_factors" : torch.tensor([1.,1.,1.,1.]),
              "allele_specific": True}
    
    # This has to be generalized and left empty 
    data_name = set(['baf', 'dr', 'dp_snp', 'vaf', 'dp'])

    def __init__(self, data_dict):
        self._params = self.params.copy()
        self._data = None
        self._name = "Clonal"
        super().__init__(data_dict, self.data_name)
        
        # internal structure is a SqueezableDict
        self._data = SqueezableDict({k:v for k,v in self._data.items()})
            
    def model(self, i = 1,  *args, **kwargs):
        n_sequences, length = 0, 0
        tot = 0

        if self._params["allele_specific"]:
            Major, minor, tot, x = self.get_Major_minor()
            #print( Major, minor, tot, x)
            
            self._params["Major"] = Major
            self._params["minor"] = minor
            probs_x = pyro.sample(
                "probs_x",
                dist.Dirichlet((1 - self._params["jumping_prob"]) * torch.eye(x.shape[0]) + self._params["jumping_prob"]).to_event(1),
            )
            
        else:
            probs_x = pyro.sample(
                "probs_x",
                dist.Dirichlet((1 - self._params["jumping_prob"]) * torch.eye(self._params["hidden_dim"]) + self._params["jumping_prob"]).to_event(1),
            )
            x = torch.arange(0, self._params["hidden_dim"]).long()
            tot = x + 1
            minor, Major = None, None
        
        if self._data["baf"] is not None:
            length, n_sequences  = self._data["baf"].shape
            baf_n_trial  = pyro.sample(
                "baf_n_trial",
                dist.Uniform(1, 10000),
            ) 
        else:
            baf_n_trial = None
            
        if self._data["dr"] is not None:
            length, n_sequences  = self._data["dr"].shape
            dr_n_trial = pyro.sample(
                "dr_n_trial",
                dist.Uniform(1, 10000),
            ) 
        else:
            dr_n_trial = None
            
        if self._data["vaf"] is not None and self._data["dp"] is not None:
            length, n_sequences  = self._data["vaf"].shape
            vaf_n_trial = pyro.sample(
                "vaf_n_trial",
                dist.Uniform(1, 10000),
            ) 
        else:
            vaf_n_trial = None
            
        if self._params["fix_purity"]:
            purity = self._params["prior_purity"]
        else:
            #purity = float(pyro.sample("purity", dist.Uniform(0.,1.)))
            purity = pyro.sample("purity", dist.Beta(4, 2))
            
                    
        if self._params["fix_ploidy"]:
            ploidy = self._params["prior_ploidy"]
        else:
            ploidy = int(pyro.sample("ploidy", dist.Poisson(2)))
            
            
        self._params["N_bins"] = length
        
        with pyro.plate("sequences", n_sequences):
            init_logits = self._params["init_probs"].log()
            trans_logits = probs_x.log()
            
            # the actual likelihood
            with ignore_jit_warnings():
                obs_dist = ClonalLikelihood(
                 x = x.unsqueeze(-1),
                 Major = Major,
                 minor = minor,
                 tot = tot,
                 baf_n_trial = baf_n_trial,
                 dr_n_trial = dr_n_trial,
                 vaf_n_trial = vaf_n_trial,
                 scaling_factors = self._params["scaling_factors"],
                 purity = purity, 
                 ploidy = ploidy,
                 batch_shape = [x.shape[0], length]
                ).to_event(1)

                hmm_dist = dist.DiscreteHMM(init_logits, trans_logits, obs_dist)
            pyro.sample("y", hmm_dist, obs=self._data)
        
    # Autoguide
    def guide(self, *args, **kwargs):
        return AutoDelta(poutine.block(self.model, hide_fn=lambda msg: msg["name"].startswith("x")))
    
    # If allele specific model map the states to the major and minor alleles
    def get_Major_minor(self):
        combinations = list(itertools.combinations_with_replacement(range( self._params["hidden_dim"]), 2))[1:]
        
        major_alleles = [max(combination) for combination in combinations]
        minor_alleles = [min(combination) for combination in combinations]

        major_allele_tensor = torch.tensor(major_alleles).long()
        minor_allele_tensor = torch.tensor(minor_alleles).long()
        x = torch.tensor(list(range(len(major_alleles)))).long()
        tot = major_allele_tensor + minor_allele_tensor

        return major_allele_tensor, minor_allele_tensor, tot,  x
    
    
    # Model 2 is the same as model 1 but with enumeration, 
    # is it used to get MAP estimates of the states, 
    # becouse Pyro is a MESSED UP library
    @infer_discrete(first_available_dim = -2, temperature=0)
    @config_enumerate
    def model_2(self,learned_params):
        
        n_sequences, length = 0, 0
        minor, Major = None, None
        x = torch.arange(1., self._params["hidden_dim"] + 1)
        tot = x
        
        if self._params["allele_specific"]:
            Major, minor, tot, x = self.get_Major_minor()
            
        probs_x = learned_params['probs_x']
        
        if self._data["baf"] is not None:
            length, n_sequences  = self._data["baf"].shape
            baf_n_trial = learned_params['baf_n_trial']
        else:
            baf_n_trial = None
            
        if self._data["dr"] is not None:
            length, n_sequences  = self._data["dr"].shape
            dr_n_trial = learned_params['dr_n_trial']
        else:
            dr_n_trial = None
            
        if self._data["vaf"] is not None:
            length, n_sequences  = self._data["vaf"].shape
            vaf_n_trial = learned_params['vaf_n_trial']
        else:
            vaf_n_trial = None
            
            
        if self._params["fix_purity"]:
            purity = self._params["prior_purity"]
        else:
            purity = learned_params['purity']
        
        
        if self._params["fix_ploidy"]:
            ploidy = self._params["prior_ploidy"]
        else:
            ploidy = learned_params['ploidy']
            
                
        with pyro.plate("sequences", n_sequences, dim = -1):
            x = [0]
            for t in pyro.markov(range(length)):
                x_new = pyro.sample(
                "x_{}".format(t),
                dist.Categorical(probs_x[x[t]])
                )
                x.append(x_new)
                
                #input_tmp = SqueezableDict({k:v[t,:] for k,v in self._data.items() if v is not None})
                
                if self._data["baf"] is not None:
                    num = (purity * minor[x_new]) +  (1 - purity)
                    den = (purity * (Major[x_new] + minor[x_new])) + (2 * (1 - purity))
                    prob = num / den
                    alpha = ((baf_n_trial-2) * prob + 1) / (1 - prob)
                    baf_lk = pyro.factor("y_baf_{}".format(t) , dist.Beta(alpha, 
                                                baf_n_trial).log_prob(self._data["baf"][t,:]))
                    # alpha = ((self._data["dp_snp"][t,:]-2) * prob + 1) / (1 - prob)
                    # baf_lk = dist.Beta(concentration1 = alpha, 
                    #                     concentration0 = self._data["dp_snp"][t,:]).log_prob(
                    #     self._data["baf"][t,:]
                    #     )
                    
                                           
                if self._data["dr"] is not None:     
                    dr = ((2 * (1-purity)) + (purity * (Major[x_new] + minor[x_new]))) / ploidy
                    dr_lk = pyro.factor("y_dr_{}".format(t) , dist.Gamma(dr * torch.sqrt(dr_n_trial) + 1, 
                                1/torch.sqrt(dr_n_trial)).log_prob(self._data["dr"][t,:]))
                
                if self._data["vaf"] is not None:
                    clonal_peaks = get_clonal_peaks(tot[x_new], Major[x_new], minor[x_new], purity)
                    tmp_vaf_lk = []
                    for j,cn in enumerate(clonal_peaks):
                        tmp_peak = 0.0
                        for i,p in enumerate(cn):
                            bin_lk = pyro.factor(f"y_vaf_{i}_{j}".format(t), dist.Binomial(self._data["dp"][t,:], 
                                    p).log_prob(self._data["vaf"][t,:].to(torch.int64)))
                    #         print(bin_lk)
                    #         tmp_peak+= (1/len(cn)) * bin_lk
                    #     tmp_vaf_lk.append(tmp_peak)
                    # vaf_lk = torch.cat(tmp_vaf_lk, dim=1)
                    
            return x
        
        
def get_clonal_peaks(tot, Major, minor, purity):
    mult = []
    for i,v in enumerate(Major):
        m = []
        if torch.equal(Major[i], minor[i]):
            m.append(Major[i][0])
        else:
            if minor[i] != 0:
                m.append(Major[i][0])
                m.append(minor[i][0])
            else:
                m.append(Major[i][0])
        if torch.equal(Major[i], torch.tensor([2])) and torch.equal(minor[i], torch.tensor([1])) == False:
            m.append(torch.tensor(1))
        mult.append(m)

    clonal_peaks = []
    for i,c in enumerate(mult):
        p = []
        for m in c:
            cp = m * purity / (tot[i] * purity + 2 * (1 - purity))
            p.append(cp)
        clonal_peaks.append(p)
        
    return clonal_peaks
