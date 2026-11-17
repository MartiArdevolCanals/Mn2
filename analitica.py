import math
import matplotlib.pyplot as plt
# PARÀMETRES I MÉS...
c_v = 3698
ro = 1081
k = 0.56
sigma = 0.472
T_in = 36.5 + 273.15  # Temperatura inicial en Kelvin
T_c = 36.5 + 273.15  # Temperatura de contorn
V = 40
h = 0.02
L = 0.02

t_0 = (c_v * ro * (L**2)) / k
T_0 = (t_0 * sigma * (V**2)) / (c_v * ro * 2 * (h**2))

T_in_norm = T_in / T_0
T_c_norm = T_c / T_0

N = 101
delta_x = 1 / (N - 1)
t_a = 0.025  # Temps adimensional d'avaluació


#DEFINICIÓ SOLUCIÓ ANALÍTICA 
# ho calculem per a un punt x, temps t 
# i ho calculem per a 50 termes de la sèrie de fourier
def s_analitica_ad(x, t, num_termes=50):
  
  sumatori = 0.0

  # Sumem només per a índexs senars (n = 1, 3, 5, ...), 
  # ja que els parells donen coeficients b_n = 0.
  for i in range(num_termes):
    n = 2 * i + 1 #índex senars

    b_n = 4.0 / (n * math.pi) # És el coeficient de Fourier b_n per a q(x) = 1

    # Part temporal (exponencial)
    terme_temporal = (1.0 - math.exp(-(n**2) * (math.pi**2) * t)) / ((n**2) * (math.pi**2))

    # Part espacial (sinus)
    terme_espacial = math.sin(n * math.pi * x)

    # Sèrie completa
    sumatori += b_n * terme_temporal * terme_espacial

  return T_c_norm + sumatori #T_c_norm actua com a base del sistema adimensional


#CÀLCUL
# Mallat espacial x de 0 a 1
x_m = [i * delta_x for i in range(N)]
x_real = [x * L * 100 for x in x_m]  # Convertit a cm (de 0 a 2 cm)

# Calculem la solució analítica per a tot el vector
T_analitica = [s_analitica_ad(x, t_a) for x in x_m]
T_real=[Ti * T_0 - 273.15 for Ti in T_analitica]

#PLOT
plt.plot(x_real, T_real)
plt.axvline(x=0.75, color="red", linestyle="--")
plt.axvline(x=1.25, color="red", linestyle="--")
#plt.axhline(y=50)
plt.show()