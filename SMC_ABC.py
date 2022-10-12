import numpy as np
from numpy.random import beta,normal

class SMC_ABC_method(object):
    def __init__(self, y, sim_params, dist_func, num_params, N, dist_final, a, c, p_acc_min):
        self.y = y
        self.sim_params = sim_params
        self.dist_func = dist_func
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

    def smc_abc_rw(self):
        part_obs = self.y
        num_drop = np.floor(self.N * self.a)
        num_keep = self.N - num_drop
        mcmc_trials = 5

        part_vals = np.zeros([self.N, self.num_params])
        part_s = np.zeros([self.N, 1])
        part_sim = np.zeros([self.N,len(part_obs)])

        for i in range(self.N):
            part_vals[i] = self.prior_sampler()
            part_sim[i] 




