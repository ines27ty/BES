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

# Tracer les grandeurs en fonction de l'angle du vilebrequin
T = [300 for x in range(0,721)]
P = [101325 for x in range(0,721)]
r = 273
m = [3.87e-5 for x in range(0,721)]

for i in range(0,len(angle)) : 
    if 0 < angle[i] < 181 :
        T[i]  = 300 
        P[i]  = 101325 
        m[i] = P[i]*V[i]/(r*T[i])*1.1
    elif 180 < angle[i] < 360 : 
        T[i] = T[i-1] *(V[i-1]/V[i])**(gamma-1)
        P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)
        m[i] = m[i-1]
    elif angle[i] == 360 : 
        T[i] = 2839.52
        P[i] = 1.13e7
        m[i] = m[i-1]
    elif 360 < angle[i] < 540 :   
        T[i] = T[i-1] *(V[i-1]/V[i])**(gamma-1) 
        P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)
        m[i] = m[i-1]
    elif angle[i] == 540 :
        T[i]  = 300
        P[i]  = 101325 
        m[i] = m[i-1]
    elif 540 < angle[i] < 720 :
        m[i] = P[i]*V[i]/(r*T[i])*1.1

   
print("Température à 0, 180, 359,360 et 540° : ", T[0], T[180], T[359], T[360], T[540])
print("Pression à 0, 180, ,359,360 et 540° : ", P[0], P[180], P[359], P[360], P[540])
print("Masse à 0, 180, 360 et 540° : ", m[0], m[180], m[360], m[540])

t = ("6 → 1 : Admission\n"
     "1 → 2 : Compression\n"
     "2 → 3 : Combustion\n"
     "3 → 4 : Détente\n"
     "4 → 5 : Échappement\n"
     "5 → 6 : Échappement")

t_angle = ("0°→180°  : Admission\n"
            "180°→360° : Compression\n"
            "360° : Combustion\n"
            "360°→540° : Détente\n"
            "540 → 720° : Échappement" )

indices = [1,2,3,4,5,6]
angles_points = [180, 359, 360, 539, 540, 720] 
V_points = [ V[180], V[359],V[360],V[539],V[540], V[720]]
P_points = [ P[180], P[359],P[360],P[539],P[540], P[720]]
m_points = [ m[180], m[359],m[360],m[539],m[540], m[720]]
T_points = [ T[180], T[359],T[360],T[539],T[540], T[720]]

indices_vol = [6,1,2,3,4,5]
angles_vol = [0, 180, 360,360, 540, 720]
points_vol = [V[0], V[180], V[360],V[360], V[540], V[720]]
pression_vol = [P[0], P[180], P[360],P[360], P[540], P[720]]
masse_vol = [m[0], m[180], m[360],m[360], m[540], m[720]]
temperature_vol = [T[0], T[180], T[360],T[360], T[540], T[720]]


plt.figure(0)
plt.plot(angle, V, label='V', color='black')
plt.scatter(angles_points, V_points, color='r', zorder=5, label='Points spécifiques',s=10)
plt.xlabel('Angle (°)')
plt.ylabel('Volume (m3)')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
for i, (x, y) in enumerate(zip(angles_vol, points_vol)):
    plt.annotate(f'{indices_vol[i]}', (x, y), textcoords="offset points", xytext=(5, 5), ha='center')
plt.title('Volume en fonction de l\'angle du vilebrequin')
plt.savefig('volume.png')
print("Volume à 0, 180, 360 et 540° : ", V[0], V[180], V[360], V[540])

# Tracer la température en fonction de l'angle du vilebrequin
plt.figure(1)
plt.plot(angle, T, label='T', color='r')
plt.scatter(angles_points, T_points, color='black', zorder=5, label='Points spécifiques',s=10)
plt.xlabel('Angle (°)')
plt.ylabel('Température (K)')
plt.grid()
plt.text(0, np.max(T)*0.7, t_angle, ha='left', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8)), fontsize=9)
rounded_positions = [round(float(pos), 1) for pos in angle]
plt.xticks(ticks=angle[::len(angle)//10 if len(angle)//10 > 0 else 1])
plt.title('Température en fonction de l\'angle du vilebrequin')
plt.savefig("T_angle_bon.png")


#Tracer le cycle thermodynamique 
plt.figure(2)
plt.plot(V, P, color='b')
plt.scatter(V_points, P_points, color='r', zorder=5, label='Points spécifiques')
plt.xlabel('Volume (m3)')
plt.ylabel('Pression (Pa)')
plt.grid()
# Annoter les points spécifiques
for i, (x, y) in enumerate(zip(V_points, P_points)):
    plt.annotate(f'{indices[i]}', (x, y), textcoords="offset points", xytext=(5, 5), ha='center')
# Ajouter une explication simple sur le graphique
plt.text(np.max(V)*0.5, np.max(P)*0.7, t, ha='left', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8)))
plt.title('Cycle thermodynamique à combustion isochore')
plt.savefig("cycle.png")



# Tracer P en fonction de l'angle du vilebrequin
plt.figure(3)
plt.plot(angle, P, label='P', color='blue')
plt.scatter(angles_points, P_points, color='black', zorder=3, label='Points spécifiques',s=10)
plt.xlabel('Angle (°)')
plt.ylabel('Pression (Pa)')
plt.grid()
plt.text(0, np.max(P)*0.7, t_angle, ha='left', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8)), fontsize=9)
plt.title('Pression en fonction de l\'angle du vilebrequin')
rounded_positions = [round(float(pos), 1) for pos in angle]
plt.xticks(ticks=angle[::len(angle)//10 if len(angle)//10 > 0 else 1])
plt.savefig("P_angle.png")


# Afficher le graphique pour la masse 
plt.figure(4)
plt.plot(angle, m, label='m', color='darkorange')
plt.scatter(angles_points, m_points, color='black', zorder=5, label='Points spécifiques',s=10)
plt.xlabel('Angle (°)')
plt.ylabel('Masse(kg)')
plt.grid()
plt.text(np.max(angle)*0.3, np.max(m)*0.2, t_angle, ha='left', bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8)), fontsize=9)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.title('Masse enfermée en fonction de l\'angle du vilebrequin')
rounded_positions = [round(float(pos), 1) for pos in angle]
plt.xticks(ticks=angle[::len(angle)//10 if len(angle)//10 > 0 else 1])
plt.savefig("m_angle.png")

#plt.show()

P_1 = 101325
Q_star = 40.3e6
T_1 = 300
cv = 1111
P_mep = P_1*(Q_star/(cv*T_1))*(1/gamma-1)*(r_c/(r_c-1))*(1-1/(r_c)**(gamma-1))
print("Pression moyenne effective : ", P_mep)


W = np.trapz(P, V) 
print("Travail total (J) : ", W)