import math
import matplotlib.pyplot as plt

#DEFINICIÓ SOLUCIÓ ANALÍTICA 
# ho calculem per a un punt x, temps t 
# fixem la temperatura del cos T_c=36,5 important per les condicions inicials i de contorn
# i ho calculem per a 50 termes de la sèrie de fourier
def solucio_analitica(x, t, T_c=0.0, num_termes=50):
  
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

  return T_c + sumatori



N = 101
dx = 1.0 / (N - 1)
t_a = 0.025  

# Mallat espacial x de 0 a 1
x_m = [i * dx for i in range(N)]

# Calculem la solució analítica per a tot el vector
T_analitica = [solucio_analitica(x, t_a) for x in x_m]

plt.plot(x_m, T_analitica)
#plt.axvline(x=0.75, color="red", linestyle="--")
#plt.axvline(x=1.25, color="red", linestyle="--")
#plt.axhline(y=50)
plt.show()