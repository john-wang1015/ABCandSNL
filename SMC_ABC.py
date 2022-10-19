'''
Sequential Monte Carlo - approximate Bayesian computation (SMC-ABC) replenishment algorithm

This Python implementation is adapt from C. C. Drovandi and A. N. Pettitt. Estimation of parameters for macroparasite
population evolution using approximate Bayesian computation.

Return:
    :param part_val:    parameter values for each particle
    :param part_sim:    summary statistics for each particle
    :param part_c:      discrepancy metric value for each particle
    :param sim:         total number of model simulations performed
    :param dist_t:      smallest discrepacy threshold reached
    :param p_acc_min:   smallest MCMC acceptance rate reached
'''

import numpy as np
from numpy.random import beta,normal
from simulator import Simulator
import scipy

class SMC_ABC_method(object):
    def __init__(self, y, sim_params, num_params, N, dist_final, a, c, p_acc_min):
        '''
        Input:
        :param y:               the observation data with the length equal to value of max_time
        :param sim_params:      a 1*3 vector contains value for page, max_time and starting Volume.
        :param num_params:      number of parameters the model have.
        :param N:               number of particle.
        :param dist_final:      target discrepancy threshold. If zero, then p_acc_min is used to determine stopping
                                criteria.
        :param a:               tuning parameter for adaptive selection of discrepancy threshold sequence.
        :param c:               tuning parameter for choosing the number of MCMC iterations in move step.
        :param p_acc_min:       minimum acceptable acceptance rate in the MCMC interations. If zero the dist_final
                                is used to determine stopping criteria.
        '''
        self.y = y
        self.sim_params = sim_params
        self.num_params = num_params
        self.N = N
        self.dist_final = dist_final
        self.a = a
        self.c = c
        self.p_acc_min = p_acc_min

    def prior_sampler(self):
        '''
            Return:
                :param sampler:     Sampler from prior distribution
        '''
        sampler = np.zeros(self.num_params)
        sampler[0] = beta(1,1)
        sampler[1] = beta(1,1e5)
        sampler[2] = np.exp(np.log(30) + normal(0,1))
        sampler[3] = np.exp(np.log(160) + normal(0,1))
        return sampler

    def dist_function(self,x,y):
        dist = np.sum(((np.log(x)) - np.log(y))**2)
        return dist

    def trans_f(self,theta):
        transf = np.zeros(self.num_params)
        transf[0] = np.log(theta[0]/(1-np.log(theta[0])))
        transf[1] = np.log(theta[1]/(1-np.log(theta[1])))
        transf[2] = np.log(theta[2]/30)
        transf[3] = np.log(theta[3]/160)
        return transf

    def pdf(self,theta_trans):
        beta_param = scipy.stats.beta.pdf(np.exp(theta_trans[0:2]/(1+ np.exp(theta_trans[0:2]))), a = [1,1], b = [1,1e5]) * np.exp(theta_trans[0:2]/(1+np.exp(theta_trans[0:2])))**2
        norm_param = scipy.stats.norm.pdf(theta_trans[2:],0,[1,1])
        prod_pdf = np.zeros(4)
        prod_pdf[0:2] = beta_param
        prod_pdf[2:] = norm_param
        theta_pdf = np.prod(prod_pdf)
        return theta_pdf

    def trans_finv(self,theta):
        finv = np.zeros(4)
        finv[0:2] = 1/(1 + np.exp(-theta[0:2]))
        finv[2] = np.exp(np.log(30) + theta[2])
        finv[3] = np.exp(np.log(160) + theta[3])
        return finv

    def smc_abc_rw(self):
        part_obs = self.y
        num_drop = np.floor(self.N * self.a)
        num_keep = int(self.N - num_drop)
        mcmc_trials = 5
        days = len(part_obs)

        part_vals = np.zeros([self.N, self.num_params])
        part_s = np.zeros([self.N, 1])
        part_sim = np.zeros([self.N,len(part_obs)])

        for i in range(self.N):
            part_vals[i] = self.prior_sampler()
            part_sim[i] = Simulator(part_vals[i][0],part_vals[i][1],part_vals[i][2].astype(np.int64),part_vals[i][3].astype(np.int64),2,part_obs[0],days).Tumourgrowth()
            part_s[i] =self.dist_function(part_obs, part_sim[i])

        sims = self.N
        dist_history = [max(part_s)]
        sims_history = [self.N]

        for i in range(self.N):
            part_vals[i] = self.trans_f(part_vals[i])

        ix = np.argsort(part_s.reshape(len(part_s)))
        part_s = np.sort(part_s).reshape(len(part_s))
        part_vals = np.asarray([part_vals[ix[i]].tolist() for i in range(len(ix))])
        part_sim = np.asarray([part_sim[ix[i]].tolist() for i in range(len(ix))])

        dist_max = part_s[self.N-1]
        dist_next = part_s[num_keep-1]
        dist_final = self.dist_final
        print(dist_max,dist_next,dist_final)
        dist_t = dist_next
        p_acc_t = 0

        while (dist_max > dist_final):
            cov_matrix = (2.38**2)*np.cov(part_vals[0:num_keep-1].T)/self.num_params

            # resample
            r = np.random.choice(num_keep, self.N - num_keep)
            part_vals[(num_keep):] = [part_vals[r[i]] for i in range(len(r))]
            part_s[(num_keep):] = [part_s[r[i]] for i in range(len(r))]
            part_sim[(num_keep):] = [part_sim[r[i]] for i in range(len(r))]

            i_acc = np.zeros(self.N-num_keep)
            sims_mcmc = np.zeros(self.N-num_keep)

            for i in range(num_keep+1,self.N):
                sum_of_dist_propr = 0
                for r in range(mcmc_trials):
                    part_vals_prop = np.random.multivariate_normal(part_vals[i],cov_matrix)
                    prior_curr = self.pdf(part_vals[i])
                    prior_prop = self.pdf(part_vals_prop)

                    #if (np.isnan(prior_prop/prior_curr) or np.random.rand() > prior_prop/prior_curr):
                    #    continue

                    #if sum_of_dist_propr <= dist_max:
                    #    continue

                    prop = self.trans_finv(part_vals_prop)
                    part_sim_prop = Simulator(prop[0],prop[1],prop[2].astype(np.int64),prop[3].astype(np.int64),2,part_obs[0],days).Tumourgrowth()
                    dist_prop = self.dist_function(part_obs,part_sim_prop)
                    sum_of_dist_propr = sum_of_dist_propr + dist_prop

                    sims_mcmc[i-num_keep] = sims_mcmc[i-num_keep]+1

                    print(dist_prop,dist_next)

                    if dist_prop <= dist_next:
                        part_vals[i] = part_vals_prop
                        part_s[i] = dist_prop
                        part_sim[i] = part_sim_prop
                        i_acc[i-num_keep] = i_acc[i-num_keep] + 1

            acc_rate = np.sum(i_acc)/(mcmc_trials*(self.N-num_keep))
            print(acc_rate)
            print(np.log(1-acc_rate))
            mcmc_iters = int(np.floor(np.log(self.c)/np.log(1-acc_rate)+1))
            print("Total number of mcmc moves for current target is {:d}, number remaining is {:d}\n".format(mcmc_iters,mcmc_iters-mcmc_trials))

            for i in range(num_keep+1,self.N):
                sum_of_dist_propr = 0
                for r in range(mcmc_iters - mcmc_trials):
                    part_vals_prop = np.random.multivariate_normal(part_vals[i],cov_matrix)
                    prior_curr = self.pdf(part_vals[i])
                    prior_prop = self.pdf(part_vals_prop)

                    if (np.isnan(prior_prop/prior_curr) or np.random.rand > prior_prop/prior_curr):
                        continue

                    if sum_of_dist_propr <= dist_max:
                        continue

                    prop = self.trans_finv(part_vals_prop)
                    part_sim_prop = Simulator(prop[0],prop[1],prop[2].astype(np.int64),prop[3].astype(np.int64),2,part_obs[0],days).Tumourgrowth()
                    dist_prop = self.dist_function(part_obs,part_sim_prop)
                    sum_of_dist_propr = sum_of_dist_propr + dist_prop

                    sims_mcmc[i-num_keep] = sims_mcmc[i-num_keep]+1

                    if dist_prop <= dist_next:
                        part_vals[i] = part_vals_prop
                        part_s[i] = dist_prop
                        part_sim[i] = part_sim_prop
                        i_acc[i-num_keep] = i_acc[i-num_keep] + 1

            num_mcmc_iters = max(0, mcmc_iters - mcmc_trials) + mcmc_trials
            p_acc = sum(i_acc)/(num_mcmc_iters*(self.N - num_keep))
            print("MCMC acceptance probability was {:d}\n".format(p_acc))

            sims = sims + sum(sims_mcmc)
            mcmc_trials = np.ceil(mcmc_iters/2)
            print("The number of unique particles is {:d}\n".format(len(np.unique(part_vals.T[0]))))

            ix = np.argsort(part_s)
            part_s = np.sort(part_s)
            part_vals = part_vals[ix]
            part_sim = part_sim[ix]

            dist_t = dist_next
            p_acc_t = p_acc

            dist_max = part_s[self.N]
            dist_next = part_s[num_keep]

            if (dist_next < dist_final):
                dist_next = dist_final

            print('The next distance is {:d} and the maximum distance is {:d} and the number to drop is {:d}\n'.format(dist_next,dist_max,num_drop))
            print('The number of sims is {:d}\n'.format(sims))

            dist_history.append(max(part_s))
            sims_history.append(sims)

            if p_acc < self.p_acc_min:
                print('Getting out as MCMC acceptance rate is below acceptable threshold\n')
                return part_vals, part_sim, part_s, sims,dist_t,p_acc_t,dist_history,sims_history

        return part_vals, part_sim, part_s, sims,dist_t,p_acc_t,dist_history,sims_history
