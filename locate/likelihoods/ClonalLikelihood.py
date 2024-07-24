import pyro
import pyro.distributions as dist
import torch
from torch.distributions import constraints
from pyro.distributions.torch_distribution import TorchDistribution
from numbers import Number



class ClonalLikelihood(TorchDistribution):
    has_rsample = False

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
                 ploidy = 2,
                 batch_shape = None,
                 validate_args = False):

        self.x = x
        self.Major = Major
        self.minor = minor
        self.tot = tot
        self.baf_n_trial = baf_n_trial
        self.dr_n_trial = dr_n_trial
        self.vaf_n_trial = vaf_n_trial
        self.scaling_factors = scaling_factors
        self.purity = purity
        self.ploidy = ploidy
        
        batch_shape = torch.Size(batch_shape)
        super(ClonalLikelihood, self).__init__(batch_shape, validate_args=validate_args)


    def log_prob(self, inp):
        dr_lk = 0
        baf_lk = 0
        vaf_lk = 0
    
        
        if self.baf_n_trial is not None:
            
            num = (self.purity * self.minor[self.x]) +  (1 - self.purity)
            den = (self.purity * (self.Major[self.x] + self.minor[self.x])) + (2 * (1 - self.purity))
            prob = num / den
            alpha = ((self.baf_n_trial-2) * prob + 1) / (1 - prob)
            baf_lk = dist.Beta(concentration1 = alpha, 
                                concentration0 = self.baf_n_trial).log_prob(
                inp["baf"]
                )
            # alpha = ((inp["dp_snp"]-2) * prob + 1) / (1 - prob)
            # baf_lk = dist.Beta(concentration1 = alpha, 
            #                     concentration0 = inp["dp_snp"]).log_prob(
            #     inp["baf"]
            #     )
            
                                    
        if self.dr_n_trial is not None:
            dr = ((2 * (1-self.purity)) + (self.purity * (self.Major[self.x] + self.minor[self.x]))) / self.ploidy
            dr_lk = dist.Gamma(dr * torch.sqrt(self.dr_n_trial) + 1, 
                                    1/torch.sqrt(self.dr_n_trial)).log_prob(
                inp["dr"]
                )
        
        if self.vaf_n_trial is not None:
            clonal_peaks = get_clonal_peaks(self.tot[self.x], self.Major[self.x], self.minor[self.x], self.purity)
            tmp_vaf_lk = []
            for cn in clonal_peaks:
                tmp_peak = 0.0
                for p in cn:
                    bin_lk = dist.Binomial(total_count = inp["dp"], 
                                                    probs = p,
                                            ).log_prob(inp["vaf"])
                    tmp_peak+= (1/len(cn)) * bin_lk
                tmp_vaf_lk.append(tmp_peak)
            vaf_lk = torch.cat(tmp_vaf_lk, dim=1)
                    
        tot_lk = self.scaling_factors[0] * baf_lk + self.scaling_factors[1] * dr_lk + self.scaling_factors[2] * vaf_lk 
        return(tot_lk)



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
