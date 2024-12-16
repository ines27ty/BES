from math import *
from math import *
import numpy as np
import matplotlib.pyplot as plt

# variables temps
dt = 0.001
tf = 10
t = 0


# conditions initiales
P0 = 101325         # Pa 
R0 = 0.001          # rayon initiale de la boule allumée
T0 = 298.15         # température initiale des gaz frais
R_p = 0.08          # à vérifier avec l'excel (rayon piston)
V = 4*np.pi*R_p**3/3  # volume du piston

rho_u = 1.2         # à calculer avec l'excel
rho_b = 2.2         # à vérifier sur l'excel
gamma = 1.3         # à vérifier sur l'excel
YF0 = 0.01          # à calculer sur l'excel    

T_gf_t = T0
T_gb_t = 2888.56    # à corriger avec la température du cerfacs
r_t = R0            # mm 
P_t = P0
QF = 1e6             # à calculer 


# Expresion de sL en fonction de phi et coefficients
phi = 1
alpha = 2.18 - 0.8*(phi-1)
beta = -0.16 + 0.22(phi-1)
B_M = 26.3                   # cm/s
B_phi = -84.7                # cm/s
phi_M = 1.13
sL0 = B_M + B_phi*(phi-phi_M)**2 #cm/s

# variables à stocker
t_tot = [T0]
P_tot = [P0]
T_gf_tot = [T0]
T_gb_tot = [2888.56]
sL_tot = [sL0]
r_tot = [R0]  


while r_t < 0.9*R_p :
    T_gb_dt = 2888.56 + 1.006*(T_gb_t - T0)     # à vérifier 
    P_dt = P_t + dt * (gamma-1)/V * QF * (4*np.pi * r_t**2 * rho_u * sL_t*YF0)
    r_dt = r_t + dt * (rho_u/rho_b) * sL0*(T_gf_t/T0)**alpha * (P_t/P0)**beta
    sL_dt = sL0 + (T_gf_t/T0)**alpha * (P_t/P0)**beta
    T_gf_dt = T_gf_t * ((R_p ** 3 - r_t ** 3) / (R_p ** 3 - r_dt ** 3)) ** (gamma - 1)


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

