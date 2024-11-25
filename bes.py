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

# Définition de la vitesse à partir de la position 
angle = np.linspace(0, 720, 721)
theta= np.radians(angle)
x = R * np.cos(np.radians(angle))
y = R * np.sin(np.radians(angle))
z = np.sqrt(x**2 + y**2)

# Calcul du volume en fonction de l'angle
term = 3.5**2 - np.sin(theta)**2

# Éviter les erreurs en imposant une valeur minimale à term (au moins 0)
V = V_c * (1 + 0.5 * (r_c - 1) * (3.5+ 1 - np.cos(theta) - np.sqrt(term)))

plt.figure(0)
plt.plot(angle, V, label='V', color='black')
plt.xlabel('Angle (°)')
plt.ylabel('Volume (m3)')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
plt.title('Volume en fonction de l\'angle du vilebrequin')
plt.savefig('volume.png')
print("Volume à 0, 180, 360 et 540° : ", V[0], V[180], V[360], V[540])


# Tracer les grandeurs en fonction de l'angle du vilebrequin
T = [300 for x in range(0,721)]
P = [101325 for x in range(0,721)]
r = 273
m = [3.87e-5 for x in range(0,721)]

for i in range(0,len(angle)) : 
    if 0 < angle[i] < 180 :
        T[i]  = 300 
        P[i]  = 101325 
        m[i] = P[i]*V[i]/(r*T[i])*1.1
    if angle[i] == 180 :
        m[i] = P[i]*V[i]/(r*T[i])*1.1
    if 180 < angle[i] < 360 : 
        T[i] = T[i-1] *(V[i-1]/V[i])**(gamma-1)
        P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)
        m[i] = m[179]
    if angle[i] == 360 : 
        T[i] = 2839.52
        P[i] = 1.13e7
        m[i] = m[359]
    if 360 < angle[i] < 540 :   
        T[i] = T[i-1] *(V[i-1]/V[i])**(gamma-1) 
        P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)
        m[i] = m[359]
    if angle[i] == 540 :
        T[i]  = 300
        P[i]  = 101325 
        m[i] = m[539]
    if 540 < angle[i] < 720 :
        m[i] = P[i]*V[i]/(r*T[i])*1.1

   
print("Température à 0, 180, 359,360 et 540° : ", T[0], T[180], T[359], T[360], T[540])
print("Pression à 0, 180, ,359,360 et 540° : ", P[0], P[180], P[359], P[360], P[540])
print("Masse à 0, 180, 360 et 540° : ", m[0], m[180], m[360], m[540])

# Afficher le graphique pour la témpérature simulée et la concentration de wo2
plt.figure(1)
plt.plot(angle, T, label='T', color='r')
plt.xlabel('Angle (°)')
plt.ylabel('Température (K)')
plt.grid()
plt.title('Température en fonction de l\'angle du vilebrequin')
rounded_positions = [round(float(pos), 1) for pos in angle]
plt.xticks(ticks=angle[::len(angle)//10 if len(angle)//10 > 0 else 1])
plt.savefig("T_angle_bon.png")

#Tracer le cycle thermodynamique 
plt.figure(2)
plt.plot(V, P, color='b')
plt.xlabel('Volume (m3)')
plt.ylabel('Pression (Pa)')
plt.grid()
plt.title('Cycle thermodynamique à combustion isochore')
plt.savefig('cycle_thermo.png')

# Tracer P en fonction de l'angle du vilebrequin
plt.figure(3)
plt.plot(angle, P, label='P', color='blue')
plt.xlabel('Angle (°)')
plt.ylabel('Pression (Pa)')
plt.grid()
plt.title('Pression en fonction de l\'angle du vilebrequin')
rounded_positions = [round(float(pos), 1) for pos in angle]
plt.xticks(ticks=angle[::len(angle)//10 if len(angle)//10 > 0 else 1])
plt.savefig("P_angle.png")


# Afficher le graphique pour la masse 
plt.figure(4)
plt.plot(angle, m, label='m', color='darkorange')
plt.xlabel('Angle (°)')
plt.ylabel('Masse(kg)')
plt.grid()
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.title('Masse enfermée en fonction de l\'angle du vilebrequin')
rounded_positions = [round(float(pos), 1) for pos in angle]
plt.xticks(ticks=angle[::len(angle)//10 if len(angle)//10 > 0 else 1])
plt.savefig("m_angle.png")

plt.show()



# Calcul du volume en fonction de l'angle du vilebrequin
# Deuxième calcul du volume
#volume = [V_d*(np.sin(np.radians(0.5*x)))+V_c for x in range(0,721)]
#volume = [V_d*z[i]*cos(np.radians(i+180))*np.pi*R**2+V_c for i in range(0,721)]
#volume = [z[i]*R**2*np.pi*np.cos(np.radians(i*0.5))+V_c for i in range(0,721)]
#volume = [V_d*np.sin(np.radians(i))+V_c for i in range(0,721)]
#V = V_c*(1+0.5*(r_c-1)*(3.5+1-np.cos(angle)-(3.5**2-(np.sin(angle))**2)**0.5))


#V = [V_c+V_d*abs(np.sin(np.radians(0.5*x))) for x in range(0,720,720//5)]
#print(V)

# tracer la fonction sinus de 0° à 720°
#V = V_d*abs(np.sin(np.radians(0.5*angle)))+V_c
#V = V_d*0.5*(1-np.cos(np.radians(angle)))+V_c
