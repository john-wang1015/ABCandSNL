import numpy as np
import matplotlib.pyplot as plt
import torch
import math
import os
from simulator import Simulator

## do parallel programming for simulated data
for i in range(10000):
    a = Simulator(0.1, 0.0001, 100, 160, 2, 200, 32).Tumourgrowth()
    if i % 1000 == 0:
        print(i)
        print(a)
    print(i)

