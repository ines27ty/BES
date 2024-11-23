from math import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter


# Tracer le volume en fonction de l'angle du vilebrequin
r_c=10
V_d = 2.78e-4
V_c = 2.78e-5
gamma = 1.28
L = 0.039
R = L/2 

# Définition de la vitesse à partir de la position 
angle = np.linspace(0, 720, 721)
x = R * np.cos(np.radians(0.5*angle))
y = R * np.sin(np.radians(0.5*angle))
z = np.sqrt(x**2 + y**2)

# Calcul du volume en fonction de l'angle du vilebrequin
#V = [z[i]*np.pi*R**2 for i in range(0,721)]


#V = [V_c+V_d*abs(np.sin(np.radians(0.5*x))) for x in range(0,720,720//5)]
#print(V)

# tracer la fonction sinus de 0° à 720°
V = V_d*abs(np.sin(np.radians(0.5*angle)))+V_c

plt.figure(0)
plt.plot(angle, V, label='V', color='r')
plt.xlabel('Angle (°)')
plt.ylabel('Volume')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
plt.title('Volume en fonction de l\'angle du vilebrequin')
plt.savefig('volume.png')


# Tracer la température en fonction de l'angle du vilebrequin
T = [300 for x in range(0,721)]
# Définition des lois de température entre 0° et 720°
for i in range(0,len(angle)) : 
    if 0 < angle[i] < 180 :
        T[i]  = 300 
    if angle[i] == 360 : 
        T[i] = 572.47
        T[i] = 2839.52
    if angle[i] == 540 :
        T[i] = 1488.02
        T[i]  = 300
    if 180 < angle[i] < 360 : 
        T[i] = T[i-1] *(V[i-1]/V[i])**(gamma-1)
    if 360 < angle[i] < 540 :   
        T[i] = T[i-1] *(V[i-1]/V[i])**(gamma-1)
   

# Afficher le graphique pour la témpérature simulée et la concentration de wo2
plt.figure(1)
fig, ax1 = plt.subplots()
ax1.plot(angle, T, label=" T", color='b')
ax1.set_xlabel("Angle (°)")
ax1.set_ylabel("Température", color='b')
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid()

# Création du deuxième axe Y
ax2 = ax1.twinx()

# Tracé de la courbe de 'Volume'
ax2.plot(angle, V, label="V", color='r')
ax2.set_ylabel("V", color='r')
ax2.tick_params(axis='y', labelcolor='r')
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

fig.tight_layout()  # Ajuste le placement pour éviter le chevauchement des labels

# Arrondir les valeurs de l'axe x à un chiffre après la virgule
rounded_positions = [round(float(pos), 1) for pos in angle]
plt.xticks(ticks=angle[::len(angle)//10 if len(angle)//10 > 0 else 1])
# Affichage de la légende combinée
fig.legend(loc='upper right', bbox_to_anchor=(0.8, 0.98))

# Sauvegarde du graphique
plt.savefig("T_V.png")


# Tracer la pression en fonction de l'angle du vilebrequin
P = [101325 for x in range(0,721)]
# Définition des lois de pression entre 0° et 720°
for i in range(0,len(angle)) : 
    if 0 < angle[i] < 180 :
        P[i]  = 101325 
    if angle[i] == 360 : 
        P[i] = 1.93e6
        P[i] = 1.13e7
    if angle[i] == 540 :
        P[i] = 5.94e5
        P[i]  = 101325 
    if 180 < angle[i] < 360 : 
        P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)        
    if 360 < angle[i] < 540 : 
        P[i] = P[i-1] *(V[i-1]/V[i])**(gamma)

    
#Tracer le cycle thermodynamique 
plt.figure(3)
plt.plot(V, P)
plt.xlabel('Volume')
plt.ylabel('Pression')
plt.grid()
plt.savefig('cycle_thermo.png')


# Afficher le graphique pour la témpérature simulée et la concentration de wo2
plt.figure(2)
fig, ax1 = plt.subplots()
ax1.plot(angle, P, label=" P", color='b')
ax1.set_xlabel("Angle (°)")
ax1.set_ylabel("Pression", color='b')
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid()

# Création du deuxième axe Y
ax2 = ax1.twinx()

# Tracé de la courbe de 'Volume'
ax2.plot(angle, V, label="V", color='r')
ax2.set_ylabel("V", color='r')
ax2.tick_params(axis='y', labelcolor='r')
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax2.grid()

fig.tight_layout()  # Ajuste le placement pour éviter le chevauchement des labels

# Arrondir les valeurs de l'axe x à un chiffre après la virgule
rounded_positions = [round(float(pos), 1) for pos in angle]
plt.xticks(ticks=angle[::len(angle)//10 if len(angle)//10 > 0 else 1])
# Affichage de la légende combinée
fig.legend(loc='upper right', bbox_to_anchor=(0.8, 0.98))

# Sauvegarde du graphique
plt.savefig("P_V.png")




# Calcul et tracé de la masse du mélange
r = 290
m = [3.85e-5 for x in range(0,721)]
# Définition des lois de température entre 0° et 720°
for i in range(0,len(angle)) : 
    if 0 < angle[i] < 720 :
        m[i] = P[i]*V[i]/(r*T[i])

# Afficher le graphique pour la masse 
plt.figure(4)
fig, ax1 = plt.subplots()
ax1.plot(angle, m, label=" m", color='b')
ax1.set_xlabel("Angle (°)")
ax1.set_ylabel("Masse (kg)", color='b')
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid()
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

# Création du deuxième axe Y
ax2 = ax1.twinx()

# Tracé de la courbe de 'Volume'
ax2.plot(angle, V, label="V", color='r')
ax2.set_ylabel("V", color='r')
ax2.tick_params(axis='y', labelcolor='r')
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

fig.tight_layout()  # Ajuste le placement pour éviter le chevauchement des labels

# Arrondir les valeurs de l'axe x à un chiffre après la virgule
rounded_positions = [round(float(pos), 1) for pos in angle]
plt.xticks(ticks=angle[::len(angle)//10 if len(angle)//10 > 0 else 1])
# Affichage de la légende combinée
fig.legend(loc='upper right', bbox_to_anchor=(0.8, 0.4))

# Sauvegarde du graphique
plt.savefig("M_V.png")

# Deuxième calcul du volume
volume = [V_d*(np.sin(np.radians(0.5*x)))+V_c for x in range(0,721)]

plt.figure(5)
plt.plot(angle, volume, label='V', color='r')
plt.xlabel('Angle (°)')
plt.ylabel('Volume')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.grid()
plt.title('Volume en fonction de l\'angle du vilebrequin')
plt.savefig('volume_cycle.png')



