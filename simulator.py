import numpy as np
from ctypes import cdll

class Simulator(object):
    def __init__(self, p0, psc, dmax, gage, page, initVolume, days):
        self.p0 = p0
        self.psc = psc
        self.dmax = dmax
        self.page = page
        self.initVolume = initVolume
        self.days = days

    def Tumourgrowth(self):
        return None
