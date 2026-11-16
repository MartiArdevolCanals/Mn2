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
N = 101
delta_x = 1/(N-1)
delta_t1 = 0.51*(delta_x**2)
delta_t2 = 0.49*(delta_x**2)
delta_t3 = 0.25*(delta_x**2)
r_1 = delta_t1/(delta_x**2)
r_2 = delta_t2/(delta_x**2)
r_3 = delta_t3/(delta_x**2)

Tn = [T_in_norm] * N
Tn1 = Tn.copy()

for i in range(1, N-1):
    Tn1[i] = Tn[i] + r_1*(Tn[i+1]-2*Tn[i]+Tn[i-1]) + delta_t1

