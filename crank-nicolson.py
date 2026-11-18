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
delta_t1 = delta_x**2
delta_t2 = 0.5*(delta_x**2)
r_1 = delta_t1/(delta_x**2)
r_2 = delta_t2/(delta_x**2)

x_norm = [i * delta_x for i in range(N)]
x = [i*L*100 for i in x_norm]
Tn = [T_in_norm] * N


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

Nt1 = round(0.025/delta_t1)
Nt2 = round(0.025)

for i in range(Nt1):
    T1=crank_nicolson(Tn, r_1, delta_t1)

for i in range(Nt2):
    T2=crank_nicolson(Tn, r_2, delta_t2)

errors1 = [abs(T1[i] - T_real[i]) for i in range(N)]
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax1.plot(x, T1)
ax1.set_xlabel("Teixit entre els electrodes (cm)")
ax1.set_ylabel("Temperatura (ºC)")
ax1.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
ax1.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
ax1.axhline(y=50, color="black", linestyle="--", label="50ºC")
ax1.legend()
ax1.set_title("Distribució de temperatura a t = 0.025 cas r=1")
ax2.plot(x, errors1, label="Error absolut", color="purple")
ax2.set_xlabel("Posició x (cm)")
ax2.set_ylabel("Error absolut (°C)")
ax2.set_title("Error respecte a la solució analítica")
plt.tight_layout()
plt.show