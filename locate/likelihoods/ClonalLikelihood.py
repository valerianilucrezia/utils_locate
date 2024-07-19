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
                 cov_atac_mu = None, 
                 cov_atac_overdispersion = None, 
                 exp_rna_mu = None,
                 exp_rna_sd = None,
                 Major = None,
                 minor = None,
                 tot = None,
                 baf_rna_number_of_trials = None,
                 baf_atac_number_of_trials = None,
                 scaling_factors = torch.tensor([1.,1.,1.,1.]),
                 purity = 1, 
                 atak_lk = "P", 
                 batch_shape = None,
                 validate_args=False):

        self.x = x
        self.cov_atac_mu = cov_atac_mu
        self.purity = purity
        self.cov_atac_overdispersion = cov_atac_overdispersion
        self.exp_rna_mu = exp_rna_mu
        self.exp_rna_sd = exp_rna_sd
        self.Major = Major
        self.minor = minor
        self.tot = tot
        self.baf_rna_number_of_trials = baf_rna_number_of_trials
        self.baf_atac_number_of_trials = baf_atac_number_of_trials
        self.scaling_factors = scaling_factors
        self.atak_lk = atak_lk

      
        batch_shape = torch.Size(batch_shape)
        
        super(MultiomeLikelihood, self).__init__(batch_shape, validate_args=validate_args)


    def log_prob(self, inp):
        
        cov_atac_lk = 0
        exp_rna_lk = 0
        baf_atac_lk = 0
        baf_rna_lk = 0
        
        
        if self.cov_atac_mu is not None:
            mean = self.purity * self.cov_atac_mu * (self.tot[self.x]) + (1 - self.purity) * self.cov_atac_mu * 2
            if self.atak_lk == "NB":
                cov_atac_lk = dist.NegativeBinomial(mean, cov_atac_overdispersion[x]).log_prob(
                    inp["cov_atac"]
                    )
            else:
                cov_atac_lk = dist.Poisson(mean).log_prob(
                                inp["cov_atac"]
                               )
        if self.exp_rna_mu is not None:
            mean = self.purity * self.exp_rna_mu[self.tot[self.x] - 1] + (1 - self.purity) * self.exp_rna_mu[1]
            
            exp_rna_lk = dist.Normal(mean, self.exp_rna_sd).log_prob(
                inp["exp_rna"]
                )
            
        if  self.baf_rna_number_of_trials is not None:
                   
            prob_tum =  (self.minor[self.x]  / (self.Major[self.x] + self.minor[self.x]))  + 1e-6
            prob = self.purity * prob_tum + 0.5 * (1 - self.purity)
            baf_rna_lk = dist.Beta(prob * self.baf_rna_number_of_trials, 
                                   (1 - prob) * self.baf_rna_number_of_trials).log_prob(
                inp["baf_rna"]
                )


        
        if  self.baf_atac_number_of_trials is not None:

            prob_tum =  (self.minor[self.x]  / (self.Major[self.x] + self.minor[self.x]))  + 1e-6
            prob = self.purity * prob_tum + 0.5 * (1 - self.purity)
            baf_atac_lk = dist.Beta(prob * self.baf_atac_number_of_trials, 
                                    (1 - prob) * self.baf_atac_number_of_trials).log_prob(
                inp["baf_atac"]
                )

        tot_lk = self.scaling_factors[0] * cov_atac_lk +  self.scaling_factors[1] * exp_rna_lk +  self.scaling_factors[2] * baf_rna_lk +  self.scaling_factors[3] * baf_atac_lk
        
        return(tot_lk)