import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import torch
import math
import os

from simulator import Simulator
from SMC_ABC import SMC_ABC_method
import SNL

y =Simulator(0.1, 0.0001, 100, 160, 2, 200, 25).Tumourgrowth()
sim_params = [2,32,200]
num_params = 4
N = 1000
dist_final = 0
a = 0.5
c = 0.01
p_acc_min = 0.2

part_vals, part_sim, part_s, sims,dist_t,p_acc_t,dist_history,sims_history = SMC_ABC_method(y,sim_params,num_params,N,dist_final,a,c,p_acc_min).smc_abc_rw()



