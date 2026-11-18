import matplotlib.pyplot as plt
from analitica import T_real

c_v = 3698
ro = 1081
k = 0.56
sigma = 0.472
T_in = 36.5+273.15
T_c = 36.5+273.15
V = 40
h = 0.02
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

def explicit(Ti, r, deltat):
    Tn1 = Ti.copy()
    for i in range(1, N-1):
        Tn1[i] = Ti[i] + r*(Ti[i+1]-2*Ti[i]+Ti[i-1]) + deltat
    return Tn1

x_norm = [i * delta_x for i in range(N)]
x = [i*L*100 for i in x_norm]

Nt1 = round(0.025/delta_t1)
T_1 = Tn.copy()
for i in range(Nt1):
    Tn1 = explicit(T_1, r_1, delta_t1)
    T_1 = Tn1.copy()
    temps1 = i*delta_t1
    
T1 = [Ti * T_0 - 273.15 for Ti in T_1]
errors1 = [abs(Tn[i] - T_real[i]) for i in range(N)]
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax1.plot(x, T1)
ax1.set_xlabel("Teixit entre els electrodes (cm)")
ax1.set_ylabel("Temperatura (ºC)")
ax1.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
ax1.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
ax1.axhline(y=50, color="black", linestyle="--", label="50ºC")
ax1.legend()
ax1.set_title("Distribució de temperatura a t = 0.025 cas r=0,51")
ax2.plot(x, errors1, label="Error absolut", color="purple")
ax2.set_xlabel("Posició x (cm)")
ax2.set_ylabel("Error absolut (°C)")
ax2.set_title("Error respecte a la solució analítica")
plt.tight_layout()
plt.show()



T_2 = Tn.copy()   
Nt2 = round(0.025/delta_t2)

for i in range(Nt2):
    Tn1 = explicit(T_2, r_2, delta_t2)
    T_2 = Tn1.copy()
    temps2 = i*delta_t2
    
T2 = [Ti * T_0 - 273.15 for Ti in T_2]
errors2 = [abs(Tn[i] - T_real[i]) for i in range(N)]
fig, (ax3, ax4) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax3.plot(x, T2)
ax3.set_xlabel("Teixit entre els electrodes (cm)")
ax3.set_ylabel("Temperatura (ºC)")
ax3.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
ax3.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
ax3.axhline(y=50, color="black", linestyle="--", label="50ºC")
ax3.legend()
ax3.set_title("Distribució de temperatura a t = 0.025 cas r=49")
ax4.plot(x, errors2, label="Error absolut", color="purple")
ax4.set_xlabel("Posició x (cm)")
ax4.set_ylabel("Error absolut (°C)")
ax4.set_title("Error respecte a la solució analítica")
plt.tight_layout()
plt.show()

T_3 = Tn.copy()
Nt3 = round(0.025/delta_t3)

for i in range(Nt3):
    Tn1 = explicit(T_3, r_3, delta_t3)
    T_3 = Tn1.copy()
    temps = i*delta_t3
    
T3 = [Ti * T_0 - 273.15 for Ti in T_3]
errors3 = [abs(Tn[i] - T_real[i]) for i in range(N)]
fig, (ax5, ax6) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax5.plot(x, T3)
ax5.set_xlabel("Teixit entre els electrodes (cm)")
ax5.set_ylabel("Temperatura (ºC)")
ax5.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
ax5.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
ax5.axhline(y=50, color="black", linestyle="--", label="50ºC")
ax5.legend()
ax5.set_title("Distribució de temperatura a t = 0.025 cas r=25")
ax6.plot(x, errors3, label="Error absolut", color="purple")
ax6.set_xlabel("Posició x (cm)")
ax6.set_ylabel("Error absolut (°C)")
ax6.set_title("Error respecte a la solució analítica")
plt.tight_layout()
plt.show()
    
plt.plot(x, T_real)
plt.show()