from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# variables temps
dt = 0.00001
tf = 10
t = 0

# conditions initiales
P0 = 101325        # pression initiale des gaz frais (Pa)
T0 = 298.15         # température initiale des gaz frais
P1 = 101325*19.3      # pression avant la combustion (Pa)
T1 = 572.47          # température avant la combustion (K)
R0 = 0.001          # rayon initiale de la boule allumée (m)
R_p = 0.05          # rayon du piston (m)
L_p = 0.03183098862 # longueur du piston (m)
#V = np.pi*R_p**2*L_p  # volume du piston (cylindre)
V = 2.82E-05        # volume des gaz dans le piston (m^3)
print("Volume du piston (m^3) = ", V)
gamma = 1.2806
YF0 = 0.0625         

# Calcul de la vitesse de flamme initiale sL0 (LIVRE)
phi = 1                             # richese
alpha = 2.18 - 0.8*(phi-1)
beta = -0.16 + 0.22*(phi-1)
B_M = 26.3*1e-2                      # m/s
B_phi = -84.7*1e-2                    # m/s
phi_M = 1.13
sL0_livre = B_M + B_phi*(phi-phi_M)**2    #m/s
#print("alpha = ", alpha)
#print("beta = ", beta)
print("sL0 (m/s) dans le livre = ", sL0_livre)            

# Nouvelle corrélation pour sL0 (THESE)
A = 55.42
B = -2.22e-14
C = -171.9
D = 74.61
E = 153.7
alpha_1 = 1.58
alpha_2 = 0.04
beta_1 = -0.203
beta_2 = -9.44e-7
sL0_ref = (A + B*(phi - phi_M) + C*(phi - phi_M)**2 + D*(phi - phi_M)**3 + E*(phi - phi_M)**4)*1e-2 # m/s
alpha_s = alpha_1 + alpha_2*(phi - phi_M)
beta_s = beta_1 + beta_2*(phi - phi_M)
sL0_these = sL0_ref * (T1/423)**alpha_s * (P1/101325)**beta_s
print("sL0_ref (m/s) = ", sL0_ref)
print("sL0 (m/s) dans la thèse= ", sL0_these)

# Calcul des densités et de r des gaz à haute température
# Gaz frais : air + c8h18 + n2 à 2000K
rho_c8h18 = 13.3
rho_air = 0.1741
rho_n2 = 0.25
y_c8h18 = 0.0625
y_air = 0.9375
y_n2 = 0.71512

rho_gf = rho_c8h18*y_c8h18 + rho_air*y_air + rho_n2*y_n2
print("rho_gf calculé = ", rho_gf)

# Gaz brûlés : h2o + co2 + n2 à 2000 K
rho_h2o = 2.09
rho_co2 = 5.12
rho_n2 = 0.259
y_h2o = 0.0885
y_co2 = 0.1923
y_n2 = 0.7192

r_mel = 8.314 / 0.03041748131
r_bru = 290.7628415
rho_gb = rho_h2o*y_h2o + rho_co2*y_co2 + rho_n2*y_n2
print("rho_gb calculé = ", rho_gb)
#rho_gf = 1.393919518
#rho_gb = 1.7803e1

# Variables en t 
sL0 = sL0_livre          # à changer si on veut celle du LIVRE ou de la THESE
sL_t = sL0               # vitesse de flamme à l'instant t
T_gf_t = T1              # température des gaz frais (avant la combustion)
T_gb_t = 2839.52         # température de fin de combustion (calculée pour isochore)
P_t = P1                 # pression à l'instant t
r_t = R0                 # m rayon initial de la boule
vb_t = 4/3*np.pi*r_t**3  # m^3 volume de la boule
QF = 44.7e6              # J/kg QLHV pouvoir calorifique du carburant
i=0                 
masse_gb = vb_t*rho_gb              # masse des gaz brûlés 
masse_gf = (V -vb_t)*rho_gf        # masse des gaz frais avec comme volume le cylindre du piston moins la boule
masse_tot = masse_gb + masse_gf                 # masse totale des gaz
print("masse totale (kg)= ", masse_tot)
print("masse gaz frais (kg)= ", masse_gf)
print("masse gaz brûlés (kg)= ", masse_gb)

# variables à stocker
t_tot = [0]
sL_tot = [sL_t]
T_gf_tot = [T1]
T_gb_tot = [2839.52]
P_tot = [P1]
r_tot = [R0]  
rho_gb_tot = [rho_gb]
rho_gf_tot = [rho_gf]
masse_gb_tot = [masse_gb]
masse_total = [masse_tot]
masse_gf_tot = [masse_gf]

while r_t < 0.9*R_p :
    T_gb_dt = 2839.52 + (1+1091.5625/1423.522842) * (T_gf_t - T0)
    P_dt = P_t + dt * (gamma-1)/V * QF * (4*np.pi * r_t**2 * rho_gf * sL_t*YF0)
    r_dt = r_t + dt * (rho_gf/rho_gb) * sL0*(T_gf_t/T0)**alpha * (P_t/P0)**beta
    #r_dt = r_t + dt * (rho_gf/rho_b) * sL_t
    sL_dt = sL0 * (T_gf_t/T0)**alpha * (P_t/P0)**beta
    T_gf_dt = T_gf_t * ((R_p ** 3 - r_t ** 3) / (R_p ** 3 - r_dt ** 3)) ** (gamma - 1)
    vb_t = 4/3*np.pi*r_t**3

    rho_gf = P_dt / (r_mel * T_gf_dt)
    rho_gb = P_dt / (r_bru * T_gb_dt)
    masse_gb = vb_t*rho_gb
    masse_gf = (V -vb_t)*rho_gf     
    #masse_gf  = rho_gf * sL_t  * 4*np.pi*r_t**2
    masse_tot = masse_gb + masse_gf
    i+=1

    # changement des valeurs de variables avant l'itération suivante
    P_t = P_dt
    T_gb_t = T_gb_dt
    r_t = r_dt
    sL_t = sL_dt
    T_gf_t = T_gf_dt
    t += dt

    # stockage des valeurs
    t_tot.append(t)
    P_tot.append(P_t)
    T_gf_tot.append(T_gf_t)
    T_gb_tot.append(T_gb_t)
    sL_tot.append(sL_t)
    r_tot.append(r_t)
    rho_gb_tot.append(rho_gb)
    rho_gf_tot.append(rho_gf)
    masse_gb_tot.append(masse_gb)
    masse_gf_tot.append(masse_gf)
    masse_total.append(masse_tot)

#print('temp gaz frais : ' + str(T_gf_tot))    
print('temps de combustion : ' + str(t))
print('itérations : ' + str(i))
print('vitesse de flamme :' + str(sL_tot[-1]))
# Tracé des courbes
plt.figure(1)
plt.loglog(t_tot, P_tot)
plt.xlabel('temps (s)')
plt.ylabel('pression (Pa)')
plt.grid()
plt.title('pression en fonction du temps')
plt.savefig('pression.png')

plt.figure(2)
plt.loglog(t_tot, T_gb_tot)
plt.xlabel('temps (s)')
plt.ylabel('température des gaz brûlés (K)')
plt.title('température des gaz brûlés en fonction du temps')
plt.grid()
plt.savefig('temp_gb.png')

plt.figure(3)
plt.loglog(t_tot, T_gf_tot)
plt.xlabel('temps (s)')
plt.ylabel('température des gaz frais (K)')
plt.title('température des gaz frais en fonction du temps')
plt.grid()
plt.savefig('temp_gf.png')

plt.figure(4)
plt.loglog(t_tot, masse_gf_tot, label='gaz frais')
plt.loglog(t_tot, masse_gb_tot, label='gaz brûlés')
plt.loglog(t_tot, masse_total, label='total')
plt.xlabel('temps (s)')
plt.ylabel('masse des gaz (kg)')
plt.title('masse des gaz en fonction du temps')
plt.legend()
plt.grid()
# Formatage en écriture scientifique pour l'axe des x
ax = plt.gca()  # Récupérer les axes actuels
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.xaxis.get_major_formatter().set_scientific(True)
ax.xaxis.get_major_formatter().set_powerlimits((0, 0))  # Force l'affichage scientifique
plt.savefig('masse_loglog.png')

plt.figure(5)
plt.plot(t_tot, masse_gf_tot, label='gaz frais')
plt.plot(t_tot, masse_gb_tot, label='gaz brûlés')
plt.plot(t_tot, masse_total, label='total')
plt.xlabel('temps (s)')
plt.ylabel('masse des gaz (kg)')
plt.title('masse des gaz en fonction du temps')
plt.legend()
plt.grid()
# Formatage en écriture scientifique pour l'axe des x
ax = plt.gca()  # Récupérer les axes actuels
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.xaxis.get_major_formatter().set_scientific(True)
ax.xaxis.get_major_formatter().set_powerlimits((0, 0))  # Force l'affichage scientifique
plt.savefig('masse.png')


plt.figure(6)
plt.loglog(t_tot, r_tot)
plt.xlabel('temps (s)')
plt.ylabel('rayon de la boule (m)')
plt.title('rayon de la boule en fonction du temps')
plt.grid()
plt.savefig('rayon.png')

plt.figure(7)
plt.loglog(t_tot, sL_tot)
plt.xlabel('temps (s)')
plt.ylabel('vitesse de flamme (m/s)')
plt.title('vitesse de flamme en fonction du temps')
plt.grid()
plt.savefig('vitesse_flamme.png')


# calcul angle nécessaire à la combustion

N_ralenti = 1200 #tr/min
N_nominal = 12000 #tr/min

angle_ralenti = 1200 / 60 * t * 360 #en °
angle_nominal = 12000 / 60 * t * 360 #en °

print('Angle nécessaire à la combustion (ralenti) : ' + str(angle_ralenti))
print('Angle nécessaire à la combustion (nominal) : ' + str(angle_nominal))

reponse = input("Afficher les graphes ? (y/n) : ").strip().lower()

if reponse == "y":
    plt.show()
else:
    print("Graphique non affiché.")