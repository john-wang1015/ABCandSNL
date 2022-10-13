import numpy as np
from numpy.random import beta,normal
from simulator import Simulator

class SMC_ABC_method(object):
    def __init__(self, y, sim_params, num_params, N, dist_final, a, c, p_acc_min):
        self.y = y
        self.sim_params = sim_params
        self.num_params = num_params
        self.N = N
        self.dist_final = dist_final
        self.a = a
        self.c = c
        self.p_acc_min

    def prior_sampler(self):
        sampler = np.zeros(self.num_params)
        sampler[0] = beta(1,1)
        sampler[1] = beta(1,1e5)
        sampler[2] = np.exp( np.log(30) + normal(0,1))
        sampler[2] = np.exp( np.log(160) + normal(0,1))
        return sampler

    def dist_function(self,x,y):
        dist = np.sum(((np.log(x)) - np.log(y))**2)
        return dist

    def trans_f(self):

        return None

    def smc_abc_rw(self):
        part_obs = self.y
        num_drop = np.floor(self.N * self.a)
        num_keep = self.N - num_drop
        mcmc_trials = 5
        days = len(part_obs)

        part_vals = np.zeros([self.N, self.num_params])
        part_s = np.zeros([self.N, 1])
        part_sim = np.zeros([self.N,len(part_obs)])

        for i in range(self.N):
            part_vals[i] = self.prior_sampler()
            part_sim[i] = Simulator.Tumourgrowth(part_vals[i][0],part_vals[i][1],part_vals[i][2],part_vals[i][3],
                                                 2,part_obs[0],days)
            part_s[i] =self.dist_function(part_obs, part_sim[i])

        sims = self.N
        dist_history = max(part_s)
        sims_history = self.N

        for i in range(self.N):
            part_vals[i]
