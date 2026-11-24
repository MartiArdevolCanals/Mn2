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
delta_t1 = delta_x**2
delta_t2 = 0.5*(delta_x**2)
r_1 = delta_t1/(delta_x**2)
r_2 = delta_t2/(delta_x**2)

#Definim el mallat espacial normalitzat i no normalitzat 
x_norm = [i * delta_x for i in range(N)]
x = [i*L*100 for i in x_norm]
#Definim la temperatura normalitzada inicial en cada punt
Tn = [T_in_norm] * N

#Definim una funció on es s'ajunta el mètode de Crank-Nicolson amb Gauss-Seidel (aquest per resoldre el sistema que ens dona Crank-Nicolson)
#Creem una llista que guarda els parametres b per cada punt (inicialment tots son 0)
#Seguidament calcula b en tots els punts 
#Després es calcula la temperatura en cada punt amb el b respectiu utilitzant la formula de Crank-Nicolson
#Aquest pas es repetix utilitzant la llista de temperatures calculada més recent, fins que la diferència entre una llista de temperatures i l'anterior és més petita que la tolerància imposada
def crank_nicolson(T, r, delta_t, tolerancia=1e-10):

    b = [0.0] * N

    for i in range(1, N-1):
        b[i] = (r/2)*T[i-1] + (1-r)*T[i] + (r/2)*T[i+1] + delta_t

    Tnova = T.copy()

    while True:
        epsilon = 0.0

        for i in range(1, N-1):

            anterior = Tnova[i]

            Tnova[i] = (b[i] + (r/2)*Tnova[i-1] + (r/2)*Tnova[i+1]) / (1+r)

            epsilon = max(epsilon, abs(Tnova[i] - anterior))

        if epsilon < tolerancia:
            break

    return Tnova

#Definim la quantitat de passos temporals
Nt1 = round(0.025/delta_t1)
Nt2 = round(0.025/delta_t2)
T1 = Tn.copy()
T2 = Tn.copy()
#Calculem la temperatura normalitzada en cada punt en el temps d'avaluació repetint la funció de Crank-Nicolson Nt vegades (ho fem pels dos casos)
for i in range(Nt1):
    T1=crank_nicolson(T1, r_1, delta_t1)
    

for i in range(Nt2):
    T2=crank_nicolson(T2, r_2, delta_t2)

#Calculem les temperatures no normalitzades en graus Celsius
T1_celsius = [temp*T_0 - 273.15 for temp in T1]
T2_celsius = [temp*T_0 - 273.15 for temp in T2]

#Calculem l'error amb la solució analítica trobada
errors1 = [abs(T1_celsius[i] - T_real[i]) for i in range(N)]

#Fem que la figura només surti quan executem aquest codi
#Fem les figures igual que pels mètodes anteriors
if __name__ == "__main__":
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    ax1.plot(x, T1_celsius)
    ax1.set_xlabel("Teixit entre els electrodes (cm)")
    ax1.set_ylabel("Temperatura (ºC)")
    ax1.tick_params(axis="x", labelbottom=True)
    ax1.axvline(x=0.75, color="red", linestyle="--", label="frontera teixit sa-malalt")
    ax1.axvline(x=1.25, color="red", linestyle="--", label="frontera teixit malalt-sa")
    ax1.axhline(y=50, color="black", linestyle="--", label="50ºC")
    ax1.legend()
    ax1.set_title(r"Distribució de temperatura $t_a$ = 0.025 cas $\gamma$=1")
    ax2.plot(x, errors1, label="Error absolut", color="purple")
    ax2.set_xlabel("Posició x (cm)")
    ax2.set_ylabel("Error absolut (°C)")
    ax2.set_title("Error respecte a la solució analítica")
    plt.tight_layout()
    plt.show()

errors2 = [abs(T2_celsius[i] - T_real[i]) for i in range(N)]
if __name__ == "__main__":
    fig, (ax3, ax4) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    ax3.plot(x, T2_celsius)
    ax3.set_xlabel("Teixit entre els electrodes (cm)")
    ax3.set_ylabel("Temperatura (ºC)")
    ax3.tick_params(axis="x", labelbottom=True)
    ax3.axvline(x=0.75, color="red", linestyle="--", label="frontera teixit sa-malalt")
    ax3.axvline(x=1.25, color="red", linestyle="--", label="frontera teixit malalt-sa")
    ax3.axhline(y=50, color="black", linestyle="--", label="50ºC")
    ax3.legend()
    ax3.set_title(r"Distribució de temperatura $t_a$ = 0.025 cas $\gamma$=0.5")
    ax4.plot(x, errors2, label="Error absolut", color="purple")
    ax4.set_xlabel("Posició x (cm)")
    ax4.set_ylabel("Error absolut (°C)")
    ax4.set_title("Error respecte a la solució analítica")
    plt.tight_layout()
    plt.show()

print("Error màxim gamma=0.5:", max(errors1), "°C")
print("Error màxim gamma=1:", max(errors2), "°C")