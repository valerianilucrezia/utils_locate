import pyro
import pyro.distributions as dist
import torch
from torch.distributions import constraints
from pyro.distributions.torch_distribution import TorchDistribution
from numbers import Number



class ClonalLikelihood(TorchDistribution):
    has_rsample = False
    #arg_constraints = {"cov_atac_mu": constraints.positive, "x" : constraints.positive, "purity": constraints.unit_interval}

    def __init__(self,
                 x = None,
                 Major = None,
                 minor = None,
                 tot = None,
                 baf_n_trial = None,
                 dr_n_trial = None,
                 vaf_n_trial = None,
                 scaling_factors = torch.tensor([1.,1.,1.,1.]),
                 purity = 1, 
                 atak_lk = "P", 
                 batch_shape = None,
                 validate_args=False):

        self.x = x
        self.Major = Major
        self.minor = minor
        self.tot = tot
        self.baf_n_trial = baf_n_trial
        self.dr_n_trial = dr_n_trial
        self.vaf_n_trial = vaf_n_trial
        self.scaling_factors = scaling_factors
        self.purity = purity
        
        batch_shape = torch.Size(batch_shape)
        
        super(ClonalLikelihood, self).__init__(batch_shape, validate_args=validate_args)


    def log_prob(self, inp):
        dr_lk = 0
        baf_lk = 0
        vaf_lk = 0
        
        if self.baf_n_trial is not None:
            prob_tum =  (self.minor[self.x]  / (self.Major[self.x] + self.minor[self.x]))  + 1e-6
            prob = self.purity * prob_tum + 0.5 * (1 - self.purity)
            baf_lk = dist.Beta(prob * self.baf_n_trial, 
                                    (1 - prob) * self.baf_n_trial).log_prob(
                inp["baf"]
                )
                                    
        if self.dr_n_trial is not None:
            dr = (self.minor[self.x] / (self.Major[self.x] + self.minor[self.x]))
            dr_lk = dist.Gamma(dr * torch.sqrt(self.dr_n_trial) + 1, 
                                    1/torch.sqrt(self.dr_n_trial)).log_prob(
                inp["dr"]
                )
        
        if self.vaf_n_trial is not None:
            pass
        
        tot_lk = self.scaling_factors[0] * baf_lk + self.scaling_factors[1] * dr_lk        
        return(tot_lk)