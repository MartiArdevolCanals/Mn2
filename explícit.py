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
Tn1 = Tn.copy()

def explicit(Ti, r, deltat):
    Tn1 = Ti.copy()
    for i in range(1, N-1):
        Tn1[i] = Ti[i] + r*(Ti[i+1]-2*Ti[i]+Ti[i-1]) + deltat
    return Tn1

x_norm = [i * delta_x for i in range(N)]
#punts_malalts = [i for i in range(N) if 0.375 <= x_norm[i] <= 0.625]
#punts_sans = [i for i in range(N) if 0.375>x_norm[i] or x_norm[i] > 0.625]

#temps = 0

#while True:

    #Tn1 = explicit(Tn, r_1, delta_t1)

    #Tmax = max(Tn1) * T_0 - 273.15
    #Tmax_sa = max(Tn1[i] for i in punts_sans) * T_0 - 273.15
    #temps = temps + delta_t1
    #Tn = Tn1.copy()

    #if temps >= 0.025:
        #break

Nt = round(0.025/delta_t2)

for i in range(Nt):
    Tn1 = explicit(Tn, r_2, delta_t2)
    Tn = Tn1.copy()
    temps = i*delta_t2
    
    

x = [i*L*100 for i in x_norm]
T = [Ti * T_0 - 273.15 for Ti in Tn]
plt.plot(x, T)
plt.axvline(x=0.75, color="red", linestyle="--")
plt.axvline(x=1.25, color="red", linestyle="--")
plt.axhline(y=50)
plt.show()
    
print(temps)