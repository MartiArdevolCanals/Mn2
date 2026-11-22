import matplotlib.pyplot as plt
import matplotlib.animation as animation
from crank_nicolson import crank_nicolson
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

punts_malalts = [i for i in range(N) if 0.375 <= x_norm[i] <= 0.625]

punts_sans = [i for i in range(N) if x_norm[i] < 0.375 or x_norm[i] > 0.625]

temps = 0

while True:

    Tn1 = crank_nicolson(Tn, r_2, delta_t2)

    Tmax = max(Tn1) * T_0 - 273.15
    Tmax_sa = max(Tn1[i] for i in punts_sans) * T_0 - 273.15

    if Tmax >= 80 or Tmax_sa >= 50:
        break

    Tn = Tn1.copy()
    temps = temps + delta_t2

T = [temp*T_0 - 273.15 for temp in Tn1]
plt.plot(x,T)
plt.axvline(x=0.75, color="red", linestyle="--", label="inici teixit malalt")
plt.axvline(x=1.25, color="red", linestyle="--", label="final teixit malalt")
plt.axhline(y=50, color="black", linestyle="--", label="50ºC")
plt.legend()
plt.show()

temps_nonorm = temps*t_0
print("Temps per a que el procediment sigui més eficaç:", temps_nonorm, "s")

# Reiniciem variables per simular i guardar els fotogrames (frames) pel GIF 
# (tot a cond inicials) temps 36,5 ; temps i pas 0
Tn = [T_in_norm] * N
temps = 0
pas = 0

#llistes per poder fer els frames totes definides i sense cap input
frames_temperatures = []
temps_registrats = []

while True:
  # Pas de Crank-Nicolson amb r_2
  Tn1 = crank_nicolson(Tn, r_2, delta_t2)

  # Comprovacions de seguretat
  Tmax = max(Tn1) * T_0 - 273.15
  Tmax_sa = max(Tn1[i] for i in punts_sans) * T_0 - 273.15

  if Tmax >= 80 or Tmax_sa >= 50:
    break

  Tn = Tn1.copy()
  temps += delta_t2

  # Desem el perfil de temperatura cada 5 passos (per evitar que el GIF vagi massa lent o pesi molt)
  #if pas % 5 == 0:

  #a cada pas ens guardem les temperatures de cada punt i el temps
  T_celsius = [temp * T_0 - 273.15 for temp in Tn1]
  frames_temperatures.append(T_celsius)
  temps_registrats.append(temps * t_0)  # Temps real en segons multiplicant per t_0
  
  pas += 1

  # Opcional: posar un límit de fotogrames si la simulació s'allarga massa
  #if len(frames_temperatures) >= 300:
  # break

# 2. Configuració de la figura de l'animació (mirar lo dels set_)

#llenç on dibuixem tot
#fig=llenç físic sencer
#ax=zona interna on es dibuixa la gràfica, amb els eixos X i Y
#figsize ample x llarg (proporcions)
fig, ax = plt.subplots(figsize=(9, 5))

#he d'acabar de mirar com funciona... li he preguntat però he hagut de marxar jeje
(line,) = ax.plot([], [], lw=2, color='blue', label='Perfil de Temperatura')

#limits eix X
ax.set_xlim(0, L * 100)
# Ajusta els límits de l'eix Y segons les temperatures que s'assoleixin (ex: de 35ºC a 85ºC)
ax.set_ylim(35,85)

#titols eixos
ax.set_xlabel('Posició $x$ a la paret cardíaca (cm)', fontsize=11)
ax.set_ylabel('Temperatura $T$ (ºC)', fontsize=11)

# Línies de referència fixes
ax.axvline(x=0.75, color='red', linestyle='--', label='Inici teixit malalt (0.75 cm)')
ax.axvline(x=1.25, color='red', linestyle='--', label='Final teixit malalt (1.25 cm)')
ax.axhline(y=50.0, color='black', linestyle=':', label='Llindar Ablació (50 ºC)')
#ens la podem petar i així el grafic es veu més gran i s'aprecia més l'evolució
ax.axhline(y=80.0, color='black', linestyle=':', label='Risc Trombosi (80 ºC)')

#llegenda
ax.legend(loc='upper right', fontsize=9)

# Text dinàmic per mostrar el temps transcorregut a cada fotograma NO CAL
t_text = ax.text(0.025,0.85,'',transform=ax.transAxes,fontsize=11,bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
#posició (0.025,0.85), comença sense res escrit (''),
# amb (transform=ax.transAxes) utilitzes coordenades relatives (0,0) inf esq  i (1,1) sup dret
#serveix pq la caixa quedi sempre en el mateix lloc sense importar l'escala dels eixos
#bbox=dict(boxstyle='round'...) defineix com és la capsa rodona...

# Funcions d'inicialització (saltable) 
def init():
 line.set_data([], [])
 t_text.set_text('')
 return line, t_text
#se suposa q abans d'executar l'evolució s'exectuta aquesta 
# és per tenir una base amb la qual treballar amb tot a 0, ben net
#ns si optimitza o q, però no es nota tant

#actualització de l'animació
def animate(i):
  line.set_data(x, frames_temperatures[i])
  t_text.set_text(f'Temps real: {temps_registrats[i]:.4f} s') #agafa 4 decimals
  return line, t_text

ani = animation.FuncAnimation(fig,animate,init_func=init,frames=len(frames_temperatures),interval=40,blit=True,)

# 3. Desar el fitxer GIF 
# (assegura't de tenir pillow instal·lat: pip install pillow)
ani.save('evolucio_ablacio.gif', writer='pillow', fps=25)

#sn el volem guardar
#plt.show()
plt.close()

