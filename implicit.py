import matplotlib.pyplot as plt
from analitica import T_real
#Mètode implícit, N=101, t_a=0.025
#\Delta t=(\Delta x)^2 o Delta t=0.5*(\Delta x)^2, \gamma=1.
#A és una matriu quadrada de 99x99
gamma1 = 0.5
gamma2 = 1
t_a = 0.025 #Falta adimensionalitzar-ho
N = 101
dx = 1/(N-1)
dt1 = 0.5*dx**2 #En aquest cas, ja que gamma=Delta t/ (Delta x)^2, i gamma = 1
num_instants1 = round(t_a/dt1)
dt2 = dx**2
num_instants2 = round(t_a/dt2)




#De la renormalització T_0=674.29 K, i com que \hat{T}=T/T_0, i inicialment T=36.5ºC=309.65 K a tots els punts, el vector temperatura comença amb tots els seus
#elements a \hat{T}=0.4592237761. Notem que la "T_vec" del codi és realment la \hat{T} adimensional.
T_ini = 0.45922
T = [T_ini]*N

#El vector b pel Gauss-Seidel d'aquest mètode d'Euler implícit té N=101 elements, dels quals el primer i últim són "T_1^{i-1}+\Delta t+gamma T_0^i", "T_{N-1}^{i-1}+\Delta t+gamma T_N^i",
#respectivament. Els elements entre mig són "T_k^{i-1}+\Delta t", per k=2,...,99. Totes aquestes "T", "t" són les versions adimensionals de les variables.
#b = [T_ini + t_a +gamma*T_ini] + [T_ini + t_a]*99 + [T_ini + t_a +gamma*T_ini]
    
def calcular_b(T_anterior, gamma, dt, T_ini):
    b = [T_anterior[1] + dt + gamma*T_ini]
    b += [T_anterior[i] + dt for i in range(2,99)]
    b += [T_anterior[99] + dt + gamma*T_ini]

    return b

A_1 = [[0.0 for j in range(99)] for i in range(99)]
A_2 = [[0.0 for j in range(99)] for i in range(99)]

for i in range(99):
    A_1[i][i] = 1 + 2*gamma1

    if i > 0:
        A_1[i][i-1] = -gamma1

    if i < 98:
        A_1[i][i+1] = -gamma1

for i in range(99):
    A_2[i][i] = 1 + 2*gamma2

    if i > 0:
        A_2[i][i-1] = -gamma2

    if i < 98:
        A_2[i][i+1] = -gamma2

def gauss_seidel(A, b_0, T_0, tol = 1e-10, iteracio_max = 1000):

    n = len(b_0)

    T = T_0.copy()
    b = b_0.copy()

    for d in range(iteracio_max):
        T_anterior = T.copy()

        # Recorrem totes les incognites
        for i in range(n):

            suma = 0.0

            for j in range(n):

                if j != i:                    
                    suma += A[i][j] * T[j]

                # Formula de Gauss-Seidel, no és la mateixa que als apunts (explicar com s'ha sintetitzat)
            T[i] = (b[i] - suma) / A[i][i] #Notem que A[i][i]=(1+2*gamma) sempre, potser es pot canviar si ha de ser més eficient

        error = max(abs(T[k]-T_anterior[k]) for k in range(n))

        if error < tol:
            return T, d+1

    print("Gauss-Seidel no ha convergit")
    return T, iteracio_max

#Falta definir el num_instants que tenim, s'ha de fer amb el t_a adimensionalitzat (t_a s'ha d'adimensionalitzar?)
T_1 = T.copy()
for i in range(num_instants1):

    #Primer s'ha de calcular b ja que es troba amb la temperatura anterior (abans d'aplicar gauss-seidel)
    b_1 = calcular_b(T_1, gamma1, dt1, T_ini)

    # Resolem amb Gauss-Seidel les 99 incògnites centrals (coneixem els extrems)
    T_centrals, iteracions = gauss_seidel(A_1, b_1, T_1[1:-1])

    #Ho apliquem a les temperatures amb els extrems fixos, això s'aplicarà al nou vector b i així successivament
    T_1 = [T_ini] + T_centrals + [T_ini]

    #print("Instant:", i+1)
    #print("Iteracions Gauss-Seidel:", iteracions)
    #print("Temperatures:", T)
T_2 = T.copy()
for i in range(num_instants2):

    #Primer s'ha de calcular b ja que es troba amb la temperatura anterior (abans d'aplicar gauss-seidel)
    b_2 = calcular_b(T_2, gamma2, dt2, T_ini)

    # Resolem amb Gauss-Seidel les 99 incògnites centrals (coneixem els extrems)
    T_centrals, iteracions = gauss_seidel(A_2, b_2, T_2[1:-1])

    #Ho apliquem a les temperatures amb els extrems fixos, això s'aplicarà al nou vector b i així successivament
    T_2 = [T_ini] + T_centrals + [T_ini]


# Representem gràficament els resultats
x = [j*dx for j in range(N)]


# --- 2. Gràfica en variables FÍSIQUES (no normalitzades) ---
x_adim = [j * dx for j in range(N)]

L_0 = 2.0        # Longitud total en cm
T_0 = 674.29     # Temperatura de referència en K

x_cm = [val * L_0 for val in x_adim]
T_celsius = [(val * T_0) - 273.15 for val in T]

