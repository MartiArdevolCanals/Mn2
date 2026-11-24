import matplotlib.pyplot as plt
from analitica import T_real

#Parametres
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

#Temperatura inicial a cada punt
Tn = [T_in_norm] * N

#Definim una funció pel mètode explícit seguint la formula que coneixem de la teoria, que es repeteix per cada punt espacial normalitzat
def explicit(Ti, r, deltat):
    Tn1 = Ti.copy()
    for i in range(1, N-1):
        Tn1[i] = Ti[i] + r*(Ti[i+1]-2*Ti[i]+Ti[i-1]) + deltat
    return Tn1

#Definim la malla espacial no normalitzada i normalitzada
x_norm = [i * delta_x for i in range(N)]
x = [i*L*100 for i in x_norm]

#Quantitat de iteracions temporals que s'han de fer, dividint el temps d'avaluació amb la delta t
Nt1 = round(0.025/delta_t1)
#Definim una funció que a repeteix el mètode d'euler explícit per cada iteració temporal  
def T_norm(T, r, deltat, Nt):
    for i in range(Nt):
        Tn1 = explicit(T, r, deltat)
        T = Tn1.copy()
    return T

#Calculem la temperatura normalitzada a cada punt en el temps d'avaluació        
T_1 = T_norm(Tn, r_1,delta_t1, Nt1)  
#Ara calculem la temperatura no normalitzada i en graus Celsius  
T1 = [Ti * T_0 - 273.15 for Ti in T_1]
#Calculem l'error a cada punt fen el valor absolut de la diferència amb el valor a cada punt trobat amb la solució analítica
errors1 = [abs(T1[i] - T_real[i]) for i in range(N)]
#Fem un figura que contingui el gràfic de la temperatura, marcant el teixit malalt i 50ºC de temperatura
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax1.plot(x, T1)
ax1.set_xlabel("Teixit entre els electrodes (cm)")
ax1.set_ylabel("Temperatura (ºC)")
ax1.tick_params(axis="x", labelbottom=True)
ax1.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
ax1.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
ax1.axhline(y=50, color="black", linestyle="--", label="50ºC")
ax1.legend()
ax1.set_title(r"Distribució de temperatura a $t_a$ = 0.025 cas $\gamma$=0,51, Euler explícit")
ax2.plot(x, errors1, label="Error absolut", color="purple")
ax2.set_xlabel("Posició x (cm)")
ax2.set_ylabel("Error absolut (°C)")
ax2.set_title("Error respecte a la solució analítica")
plt.tight_layout()
plt.show()

#Repetim per els altres dos casos
Nt2 = round(0.025/delta_t2)
T_2 = T_norm(Tn, r_2, delta_t2, Nt2)
T2 = [Ti * T_0 - 273.15 for Ti in T_2]
errors2 = [abs(T2[i] - T_real[i]) for i in range(N)]
fig, (ax3, ax4) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax3.plot(x, T2)
ax3.set_xlabel("Teixit entre els electrodes (cm)")
ax3.set_ylabel("Temperatura (ºC)")
ax3.tick_params(axis="x", labelbottom=True)
ax3.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
ax3.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
ax3.axhline(y=50, color="black", linestyle="--", label="50ºC")
ax3.legend()
ax3.set_title(r"Distribució de temperatura a $t_a$ = 0.025 cas $\gamma$=0.49, Euler explícit")
ax4.plot(x, errors2, label="Error absolut", color="purple")
ax4.set_xlabel("Posició x (cm)")
ax4.set_ylabel("Error absolut (°C)")
ax4.set_title("Error respecte a la solució analítica")
plt.tight_layout()
plt.show()


Nt3 = round(0.025/delta_t3)
T_3 = T_norm(Tn, r_3, delta_t3, Nt3)
T3 = [Ti * T_0 - 273.15 for Ti in T_3]
errors3 = [abs(T3[i] - T_real[i]) for i in range(N)]
fig, (ax5, ax6) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax5.plot(x, T3)
ax5.set_xlabel("Teixit entre els electrodes (cm)")
ax5.set_ylabel("Temperatura (ºC)")
ax5.tick_params(axis="x", labelbottom=True)
ax5.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
ax5.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
ax5.axhline(y=50, color="black", linestyle="--", label="50ºC")
ax5.legend()
ax5.set_title(r"Distribució de temperatura $t_a$ = 0.025 cas $\gamma$=0.25, Euler explícit")
ax6.plot(x, errors3, label="Error absolut", color="purple")
ax6.set_xlabel("Posició x (cm)")
ax6.set_ylabel("Error absolut (°C)")
ax6.set_title("Error respecte a la solució analítica")
plt.tight_layout()
plt.show()
    

#Error màxim per cada cas
print("Error màxim cas 1:", max(errors1))
print("Error màxim cas 2:", max(errors2))
print("Error màxim cas 3:", max(errors3))