from math import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
from scipy.ndimage import gaussian_filter1d

# Tracer le volume en fonction de l'angle du vilebrequin
r_c=10
V_d = 2.78e-4
V_c = 2.78e-5
gamma = 1.28
L = 0.039
R = L/2 

# Définition de la vitesse à partir de l'angle du vilebrequin 
angle = np.linspace(0, 720, 721)
theta= np.radians(angle)

term = 3.5**2 - np.sin(theta)**2
V = V_c * (1 + 0.5 * (r_c - 1) * (3.5+ 1 - np.cos(theta) - np.sqrt(term)))


# Tracer les valeurs pour tout le cycle en fonction de l'angle du vilebrequin
r = 273
rho = 1.37e1

P = [101325 for x in range(0,721)]
T = [300 for x in range(0,721)]
m = [3.87e-5 for x in range(0,721)]
cv = 1111.57
E = -2800000
volume = [0 for x in range(0,721)]

for i in range(0,len(angle)) : 
    if 0 < angle[i] < 181 :
        P[i]  = 101325 
        T[i] = 300
        m[i] = P[i]*V[i]/(r*T[i])*1.1
    elif 180 < angle[i] < 360 :
        P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)
        T[i] = T[i-1] * (V[i-1]/V[i])**(gamma-1)
        m[i] = P[i]*V[i]/(r*T[i])*1.1
    elif angle[i] == 360 :
        P[i] = 7.5e6
        m[i] = m[i-1]
        T[i] = P[i]*V[i]/(290*m[i])
    elif 360 < angle[i] < 541 :
        m[i] = m[360]
        P[i] = P[360]
        m_brulee = (T[i-1] - T[i])* m[i]*cv/E 
        m_restante = m[i] - m_brulee
        #print("Masse brûlée à 360° : ", m_brulee)
        #print("Masse restante à 360° : ", m_restante)
        T[i] = T[i-1] - m_restante*E/(m[i]*cv)
        volume[i] = m[i]*r*T[i]/P[i]
        T[i] = P[i]*V[i]/(r*m[i])
        if V[i] > 4.24e-5 :
            P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)
            T[i] = T[i-1] *(V[i-1]/V[i])**(gamma-1)
    elif angle[i] == 540 :
        m[i]= m[i-1]
    elif 540 < angle[i] < 720 :
        m[i] = P[i]*V[i]/(r*T[i])*1.1


t = ("6 → 1 : Admission\n"
     "1 → 2 : Compression\n"
     "2 → 3 : Combustion isoV\n"
     "3' → 3 : Combustion isoP\n"
     "3 → 4 : Détente\n"
     "4 → 5 : Échappement\n"
     "5 → 6 : Échappement")

t_angle = ("0°→180°  : Admission\n"
            "180°→360° : Compression\n"
            "360°→385° : Comb. isoV et isoP\n"
            "385°→540° : Détente\n"
            "540 → 720° : Échappement" )

indices = [1,2,3,"3'",4,5,6]
V_points = [ V[180], V[359],V[360],V[384],V[539],V[541], V[720]]
P_points = [ P[180], P[359],P[360],P[384],P[539],P[541], P[720]]

indices_angles = [6,1,2,3,"3'",4,5]
points_angles = [0,180, 359,360, 384, 540, 720]
volume_angles = [ V[0],V[180], V[359],V[360],V[384],V[540], V[720]]
pression_angles = [ P[0],P[180], P[359],P[360],P[384],P[540], P[720]]
masse_angles = [ m[0],m[180], m[359],m[360],m[384],m[540], m[720]]
temperature_angles = [ T[0],T[180], T[359],T[360],T[384],T[540], T[720]]

plt.figure(0)
plt.plot(angle, V, label='V', color='black')
plt.scatter(points_angles, volume_angles, label='V', color='red')
plt.xlabel('Angle (°)')
plt.ylabel('Volume (m3)')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
for i, (x, y) in enumerate(zip(points_angles, volume_angles)):
    plt.annotate(f'{indices_angles[i]}', (x, y), textcoords="offset points", xytext=(5, 5), ha='center')
plt.title('Volume en fonction de l\'angle du vilebrequin')
plt.savefig('isob_volume.png')


plt.figure(1)
plt.plot(angle, P, label='P', color='b')
plt.scatter(points_angles, pression_angles, label='P', color='red')
plt.xlabel('Angle (°)')
plt.ylabel('Pression (Pa)')
plt.text(0, np.max(P)*0.7, t_angle, ha='left', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8)), fontsize=9)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
for i, (x, y) in enumerate(zip(points_angles, pression_angles)):
    plt.annotate(f'{indices_angles[i]}', (x, y), textcoords="offset points", xytext=(5, 5), ha='center')
plt.title('Pression en fonction de l\'angle du vilebrequin')
plt.savefig('isob_pressure.png')

print("Pression à 0, 180, 359, 360 et 540° : ", P[0], P[180], P[359], P[360], P[540])
print("Volume à 0, 180, 359, 360 et 540° : ", V[0], V[180], V[359], V[360], V[540])
print("Température à 0, 180, 359, 360 et 540° : ", T[0], T[180], T[359], T[360], T[540])

plt.figure(2)
plt.plot(angle, T, label='T', color='red')
plt.scatter(points_angles, temperature_angles, label='T', color='red')
plt.xlabel('Angle (°)')
plt.ylabel('Température (K)')
plt.grid()
plt.text(0, np.max(T)*0.8, t_angle, ha='left', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8)), fontsize=9)
for i, (x, y) in enumerate(zip(points_angles, temperature_angles)):
    plt.annotate(f'{indices_angles[i]}', (x, y), textcoords="offset points", xytext=(5, 5), ha='center')
plt.title('Température en fonction de l\'angle du vilebrequin')
plt.savefig('isob_temperature.png')
print("Température à 0, 180, 360 et 540° : ", T[0], T[180], T[360], T[540])

plt.figure(3)
plt.plot(angle, m, label='m', color='darkorange')
plt.scatter(points_angles, masse_angles, label='m', color='red')
plt.xlabel('Angle (°)')
plt.ylabel('Masse (kg)')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
plt.text(np.max(angle)*0.3, np.max(m)*0.2, t_angle, ha='left', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8)), fontsize=9)
for i, (x, y) in enumerate(zip(points_angles, masse_angles)):
    plt.annotate(f'{indices_angles[i]}', (x, y), textcoords="offset points", xytext=(5, 5), ha='center')
plt.title('Masse en fonction de l\'angle du vilebrequin')
plt.savefig('isob_mass.png')
print("Masse à 0, 180, 360 et 540° : ", m[0], m[180], m[360], m[540])


# Tracer le cycle thermodynamique
plt.figure(4)
plt.plot(V, P, color='b')
plt.scatter(volume_angles, pression_angles, color='red')
plt.xlabel('Volume (m3)')
plt.ylabel('Pression (Pa)')
plt.grid()
plt.text(np.max(V)*0.5, np.max(P)*0.7, t, ha='left', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8)))
for i, (x, y) in enumerate(zip(V_points, P_points)):
    plt.annotate(f'{indices[i]}', (x, y), textcoords="offset points", xytext=(5, 5), ha='center')
plt.title('Cycle thermodynamique à combustion isochore puis isobare')
plt.savefig('isob_cycle.png')

plt.figure(5)
plt.plot(angle,V, color='red')
plt.plot(angle,volume, color='red')
plt.xlabel('Angle (°)')
plt.ylabel('Volume (m3)')
plt.grid()


# Calcul énergétique : 

#plt.show()
alpha = P[360] / P[359]
beta =  V[384] / V[359]
eta_cycle = 1 - (1 / r_c**(gamma - 1)) * ((alpha * beta**gamma - 1) / (alpha * gamma * (beta - 1) + alpha - 1))
print("Alpha : ", alpha)
print("Beta : ", beta)
print("Efficacité du cycle : ", eta_cycle)

P_1 = 101325
T_1 = 300
cv = 1111
c_p= 1432
Q_star = cv*(T[360] - T[359]) + c_p*(T[384] - T[360])
print("Q* (J/kg) : ", Q_star)

P_mep = P_1 * (Q_star / (cv * T_1 * (gamma - 1))) * ((r_c) / (r_c - 1)) * eta_cycle
P_3 = P[360]

P_mep_bis = P_3 * (Q_star / (cv * T_1 * (gamma - 1))) * ((r_c) / (r_c - 1)) * eta_cycle / (alpha * r_c ** gamma)
print("Pmep (bar): ", P_mep/1e5)
print("Pmep bis (bar): ", P_mep_bis/1e5)
W = np.trapz(P, V) 
print("Travail total (J) : ", W)
