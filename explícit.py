import matplotlib.pyplot as plt

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

for i in range(Nt1):
    Tn1 = explicit(Tn, r_1, delta_t1)
    Tn = Tn1.copy()
    temps1 = i*delta_t1
    
T1 = [Ti * T_0 - 273.15 for Ti in Tn]
plt.plot(x, T1)
plt.xlabel("Teixit entre els electrodes (cm)")
plt.ylabel("Temperatura (ºC)")
plt.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
plt.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
plt.axhline(y=50, color="black", linestyle="--", label="50ºC")
plt.legend()
plt.show()

Tn = [T_in_norm] * N    
Nt2 = round(0.025/delta_t2)

for i in range(Nt2):
    Tn1 = explicit(Tn, r_2, delta_t2)
    Tn = Tn1.copy()
    temps2 = i*delta_t2
    
T2 = [Ti * T_0 - 273.15 for Ti in Tn]
plt.plot(x, T2)
plt.xlabel("Teixit entre els electrodes (cm)")
plt.ylabel("Temperatura (ºC)")
plt.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
plt.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
plt.axhline(y=50, color="black", linestyle="--", label="50ºC")
plt.legend()
plt.show()

Tn = [T_in_norm] * N 
Nt3 = round(0.025/delta_t3)

for i in range(Nt3):
    Tn1 = explicit(Tn, r_3, delta_t3)
    Tn = Tn1.copy()
    temps = i*delta_t3
    
T3 = [Ti * T_0 - 273.15 for Ti in Tn]
plt.plot(x, T3)
plt.xlabel("Teixit entre els electrodes (cm)")
plt.ylabel("Temperatura (ºC)")
plt.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
plt.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
plt.axhline(y=50, color="black", linestyle="--", label="50ºC")
plt.legend()
plt.show()
    
