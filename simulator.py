import numpy as np
from ctypes import cdll, c_double,c_int
import ctypes

class Simulator(object):
    def __init__(self, p0, psc, dmax, gage, page, initVolume, days):
        self.p0 = p0
        self.psc = psc
        self.dmax = dmax
        self.gage = gage
        self.page = page
        self.initVolume = initVolume
        self.days = days

    def Tumourgrowth(self):
        lib = cdll.LoadLibrary('./VCBM/Model.so')
        lib.TumourGrowthData.argtypes = [ctypes.POINTER(ctypes.c_double),c_double,c_double,c_int,c_int,c_int,c_double,c_int]
        input = ctypes.ARRAY(ctypes.c_double,32)()
        pf = lib.TumourGrowthData(input, self.p0, self.psc, self.dmax, self.gage, self.page, self.initVolume, self.days)
        data = np.zeros(self.days)

        for i in range(self.days):
            data[i] = input[i]

        return data

if __name__ == "__main__":
    '''
        A test for simulator
    '''
    simulator = Simulator(0.1, 0.0001, 100, 160, 2, 200, 32)
    print(simulator.Tumourgrowth())
