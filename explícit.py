import numpy as np
import matplotlib.pyplot as plt

c_v = 3698
ro = 1081
k = 0.56
sigma = 0.472
T_in = 36.5+273.15
T_c = 36.5+273.15
V = 40
h = 0.005
L = 0.02
t_0 = (c_v*ro*(L**2))/k
T_0 = (t_0*sigma*(V**2))/(c_v*ro*2*(h**2))
T_in_norm = T_in/T_0
T_c_norm = T_c/T_0