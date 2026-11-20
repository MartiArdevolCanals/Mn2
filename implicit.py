import matplotlib.pyplot as plt
from analitica import T_real
#Mètode implícit, N=101, t_a=0.025
#\Delta t=(\Delta x)^2, \gamma=1.
#A és una matriu quadrada de 99x99
gamma = 1
t_a = 0.025 #Falta adimensionalitzar-ho
N = 101
dx = 1/(N-1)
dt = dx**2 #En aquest cas, ja que gamma=Delta t/ (Delta x)^2, i gamma = 1
num_instants = round(t_a/dt)



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

A = [[0.0 for j in range(99)] for i in range(99)]

for i in range(99):
    A[i][i] = 1 + 2*gamma

    if i > 0:
        A[i][i-1] = -gamma

    if i < 98:
        A[i][i+1] = -gamma

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

for i in range(num_instants):

    #Primer s'ha de calcular b ja que es troba amb la temperatura anterior (abans d'aplicar gauss-seidel)
    b = calcular_b(T, gamma, dt, T_ini)

    # Resolem amb Gauss-Seidel les 99 incògnites centrals (coneixem els extrems)
    T_centrals, iteracions = gauss_seidel(A, b, T[1:-1])

    #Ho apliquem a les temperatures amb els extrems fixos, això s'aplicarà al nou vector b i així successivament
    T = [T_ini] + T_centrals + [T_ini]

    #print("Instant:", i+1)
    #print("Iteracions Gauss-Seidel:", iteracions)
    #print("Temperatures:", T)


# Representem gràficament els resultats
x = [j*dx for j in range(N)]

plt.plot(x, T)
plt.xlabel("Posicio x adimensional")
plt.ylabel("Temperatura adimensional")
plt.title("Euler implicit")
plt.grid()
plt.show()

# --- 2. Gràfica en variables FÍSIQUES (no normalitzades) ---
x_adim = [j * dx for j in range(N)]

L_0 = 2.0        # Longitud total en cm
T_0 = 674.29     # Temperatura de referència en K

x_cm = [val * L_0 for val in x_adim]
T_celsius = [(val * T_0) - 273.15 for val in T]

plt.subplot(1, 2, 2)
plt.plot(x_cm, T_celsius, color='red')
plt.xlabel("Posició x (cm)")
plt.ylabel("Temperatura (°C)")
plt.title("Variables físiques (reals)")
plt.grid(True)

plt.tight_layout()
plt.show()