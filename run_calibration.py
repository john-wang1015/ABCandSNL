import numpy as np
import matplotlib.pyplot as plt
from simulator import Simulator
from SMC_ABC import SMC_ABC_method
import SNL
import torch

# SMC-ABC
y = Simulator(0.1, 0.0001, 100, 160, 2, 200, 5).Tumourgrowth()
sim_params = [2,32,200]
num_params = 4
N = 1000
dist_final = 0
a = 0.5
c = 0.01
p_acc_min = 0.01

part_vals, part_sim, part_s, sims,dist_t,p_acc_t,dist_history,sims_history = SMC_ABC_method(y,sim_params,num_params,N,dist_final,a,c,p_acc_min).Sampler()

#  SNL
N_iter = 2
N_sim = 10000
D_sams = np.zeros([N_iter*N_sim, num_params])

for i in range(N_sim):
    D_sams[i] = SNL.CustomPriorDist(30*torch.ones(1),160*torch.ones(1),torch.zeros(1),torch.ones(1),torch.ones(1)*1e4).sample().numpy()

D_sams[N_sim:N_iter*N_sim] = SNL.SNL(D_sams[0:N_sim],N_sim,y).infer_and_sanpler()
part_sim_snl = np.zeros([N_sim, len(y)])

for i in range(N_sim):
    part_sim_snl[i] = Simulator(float(D_sams[N_sim+i][0]), float(D_sams[N_sim+i][1]), int(D_sams[N_sim+i][2]), int(D_sams[N_sim+i][3]), 2, float(y[0]), int(len(y))).Tumourgrowth()

