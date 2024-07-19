""" Utils class

A set of utils function to run automatically an enetire inference cycle, plotting and saving results.

"""

import matplotlib.pyplot as plt
import pandas as pd
import os
from locate.inference import *
import numpy as np
from pyro.optim import ClippedAdam
from pyro.infer import SVI, TraceGraph_ELBO
import torch
from pandas.core.common import flatten

def plot_loss(loss, save = False, output = "run1"):
    plt.plot(loss)
    plt.title("ELBO")
    plt.xlabel("step")
    plt.ylabel("loss")
    if(save):
        plt.savefig(output + "_ELBO.png")

def collect_params(pars):
    pars = list(flatten([value.detach().tolist() for key, value in pars.items()]))
    return(np.array(pars))

def retrieve_params():
    param_names = pyro.get_param_store()
    res = {nms: pyro.param(nms) for nms in param_names}
    return res
