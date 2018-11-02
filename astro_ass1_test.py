import numpy as np
import math
import time as tm
import matplotlib.pyplot as plt

N = 6  # at the moment only even values of N allowed due to initialisation
w = 4
dt = 0.1
dx = 0.2
rho = 1.05
v = 2.05
P = 3.05
gamma = 1.4

for i in range(N+2):
    print(i)
