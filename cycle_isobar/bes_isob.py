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

plt.figure(0)
plt.plot(angle, V, label='V', color='black')
plt.xlabel('Angle (°)')
plt.ylabel('Volume (m3)')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
plt.title('Volume en fonction de l\'angle du vilebrequin')
plt.savefig('isob_volume.png')
print("Volume à 0, 180, 360 et 540° : ", V[0], V[180], V[360], V[540])

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
    if 180 < angle[i] < 360 :
        P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)
        T[i] = T[i-1] * (V[i-1]/V[i])**(gamma-1)
        m[i] = P[i]*V[i]/(r*T[i])*1.1
    if angle[i] == 360 :
        P[i] = 7.5e6
        m[i] = m[i-1]
        T[i] = P[i]*V[i]/(290*m[i])
       
    if 360 < angle[i] < 541 :
        m[i] = m[360]
        P[i] = P[360]
        m_brulee = (T[i-1] - T[i])* m[i]*cv/E 
        m_restante = m[i] - m_brulee
        print("Masse brûlée à 360° : ", m_brulee)
        print("Masse restante à 360° : ", m_restante)
        T[i] = T[i-1] - m_restante*E/(m[i]*cv)
        volume[i] = m[i]*r*T[i]/P[i]
        T[i] = P[i]*V[i]/(r*m[i])
        if V[i] > 4.24e-5 :
            P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)
            T[i] = T[i-1] *(V[i-1]/V[i])**(gamma-1)
    if angle[i] == 540 :
        m[i]= m[i-1]
    if 540 < angle[i] < 720 :
        m[i] = P[i]*V[i]/(r*T[i])*1.1
plt.figure(1)
plt.plot(angle, P, label='P', color='b')
plt.xlabel('Angle (°)')
plt.ylabel('Pression (Pa)')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
plt.title('Pression en fonction de l\'angle du vilebrequin')
plt.savefig('isob_pressure.png')
print("Pression à 0, 180, 360 et 540° : ", P[0], P[180], P[360], P[540])
print("Température à 359,360,361, 540° : ", T[359], T[360], T[361], T[540]) 

plt.figure(2)
plt.plot(angle, T, label='T', color='red')
plt.xlabel('Angle (°)')
plt.ylabel('Température (K)')
plt.grid()
plt.title('Température en fonction de l\'angle du vilebrequin')
plt.savefig('isob_temperature.png')
print("Température à 0, 180, 360 et 540° : ", T[0], T[180], T[360], T[540])

plt.figure(3)
plt.plot(angle, m, label='m', color='darkorange')
plt.xlabel('Angle (°)')
plt.ylabel('Masse (kg)')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
plt.title('Masse en fonction de l\'angle du vilebrequin')
plt.savefig('isob_mass.png')
print("Masse à 0, 180, 360 et 540° : ", m[0], m[180], m[360], m[540])


# Tracer le cycle thermodynamique
plt.figure(4)
plt.plot(V, P, color='b')
plt.xlabel('Volume (m3)')
plt.ylabel('Pression (Pa)')
plt.grid()
plt.title('Cycle thermodynamique à combustion isochore puis isobare')
plt.savefig('isob_cycle.png')

plt.figure(5)
plt.plot(angle,V, color='red')
plt.plot(angle,volume, color='red')
plt.xlabel('Angle (°)')
plt.ylabel('Volume (m3)')
plt.grid()