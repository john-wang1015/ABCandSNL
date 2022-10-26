import numpy as np
from simulator import Simulator
import torch
_ = torch.manual_seed(10)
import math

import sbi
from sbi import utils as utils
from sbi.inference import SNLE, prepare_for_sbi, simulate_for_sbi
from sbi.utils import process_prior

class CustomPriorDist(object):
    def __init__(self, loc1, loc2, mean,sigma, psc, return_numpy: bool = False):
        self.loc1 = loc1
        self.loc2 = loc2
        self.psc = psc
        self.dist1 = torch.distributions.normal.Normal(mean, sigma)
        self.dist2 = utils.BoxUniform(mean,sigma)
        self.dist3 = torch.distributions.beta.Beta(sigma,psc)
        self.return_numpy = return_numpy

    def sample(self, sample_shape=torch.Size([])):
        if len(sample_shape) == 1:
            length = sample_shape[0]
            samples = torch.ones(length,4)
            temp_p0 = self.dist2.sample(sample_shape)
            temp_psc = self.dist3.sample(sample_shape)
            temp_dmax = torch.exp(torch.log(self.loc1*torch.ones(sample_shape)) + self.dist1.sample(sample_shape))
            temp_gage = torch.exp(torch.log(self.loc2*torch.ones(sample_shape)) + self.dist1.sample(sample_shape))
            samples[:,0] = temp_p0[:,0]
            samples[:,1] = temp_psc[:,0]
            samples[:,2] = temp_dmax[:,0]
            samples[:,3] = temp_gage[:,0]
            return samples.numpy() if self.return_numpy else samples
        else:
            samples = torch.ones(1,4)
            temp_p0 = self.dist2.sample(sample_shape)
            temp_psc = self.dist3.sample(sample_shape)
            temp_dmax = torch.exp(torch.log(self.loc1*torch.ones(sample_shape)) + self.dist1.sample(sample_shape))
            temp_gage = torch.exp(torch.log(self.loc2*torch.ones(sample_shape)) + self.dist1.sample(sample_shape))
            samples[0,0] = temp_p0[0]
            samples[0,1] = temp_psc[0]
            samples[0,2] = torch.abs(temp_dmax[0])
            samples[0,3] = temp_gage[0]
            return samples.numpy() if self.return_numpy else samples

    def log_prob(self, values):
        log_probs = torch.ones((values.size()[0],))
        length = values.size()[0]
        if self.return_numpy:
            values = torch.as_tensor(values)

        for i in range(values.size()[0]):
            temp = torch.ones(4)
            temp[0] = self.dist2.log_prob(values[i][0])
            temp[1] = self.dist3.log_prob(torch.abs(values[i][1]))
            temp[2] = torch.log(torch.exp(torch.log(self.loc1*torch.ones(1)) + 1/(torch.sqrt(2*math.pi*torch.ones(1)))*torch.exp(-0.5 * (values[i][2]-0)**2)))
            temp[3] = torch.log(torch.exp(torch.log(self.loc2*torch.ones(1)) + 1/(torch.sqrt(2*math.pi*torch.ones(1)))*torch.exp(-0.5 * (values[i][3]-0)**2)))
            log_probs[i] =  torch.sum(temp)

        return log_probs.numpy() if self.return_numpy else log_probs

class SNL(object):
    def __init__(self, theta,n_samples,observation):
        self.theta = theta
        self.n = n_samples
        self.observation = observation
        self.simulation = np.zeros(self.n, len(self.observation))

        for i in range(self.n):
            self.simulation[i] = Simulator(self.theta[0],self.theta[1],self.theta[2].astype(np.int64),self.theta[3].astype(np.int64),2,self.observation[0],len(self.observation)).Tumourgrowth()

        priorDist = CustomPriorDist(30*torch.ones(1),160*torch.ones(1),torch.zeros(1),torch.ones(1),torch.ones(1)*1e4)
        self.prior = process_prior(priorDist)[0]

    def processing(self):
        theta = torch.from_numpy(self.theta).to(torch.float32)
        x = torch.from_numpy(self.simulation).to(torch.float32)
        x_0 = torch.from_numpy(self.observation).to(torch.float32)
        return theta,x,x_0

    def infer_and_sanpler(self):
        inference = sbi.inference.SNLE(prior=self.prior, density_estimator='nsf')
        theta,x,x_0 = self.processing()
        density_estimator = inference.append_simulations(theta, x).train()
        posterior = inference.build_posterior(density_estimator)
        posterior_samples = posterior.sample((self.n,),x = x_0)
        return posterior_samples
