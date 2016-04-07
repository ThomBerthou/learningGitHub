# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 14:09:54 2014

@author: Thomas B
Cette fonction permet de calculer le flux solaire direct, diffus et réfléchi arrivant sur les parois du bâtiment. 
Dans la V0.3 il n’y a que 4 orientations possibles pour les façades : nord, sud, est, ouest. Les algorithmes proviennent de la thèse de  David Da Silva (2011).

GMT=nombre d'heure de décalage avec le méridien de greenwich (fixé à 1)
rho=albédo de surface (fixé à 0.2)
beta=inclinaison de la surface (ficcée à 90°)
gamma=azimuth de la surface(en degré,sud=0°est=-90°ouest=90°nord=-180 ou +180°
"""
import numpy as np
import transformer as tr
def solar_flux(city,longitude, latitude, I_dif, I_dir, delta):
    
    NB_dwe = len(city)
    pi=3.14159
    GMT=1.
    nb_hour=np.linspace(1,8760,8760)
    deg2rad = pi/180. #conversion degré en radian
    rad2deg = 180./pi #conversion radian en degré
    heure_23 = np.tile(np.linspace(0,23,24),365)
    
    rho = 0.2
    beta = 90
    
# Declinaison  solaire, n est le numero de jour
    n = np.around(nb_hour / 24 + 0.5)
    n=np.r_[1,n]
    n=n[0:8760]
    declinaison=deg2rad * 23.45 * np.sin((deg2rad) * 360 * (284 + n) / 365.)
    
#Décalage à l'heure civile
    B = ((n-1) * 360 / 365) * deg2rad
    E = 229.2 * (0.000075 + 0.001868 * np.cos(B) - 0.032077 * np.sin(B) - 0.014615
    * np.cos(2*B) - 0.04089 * np.sin(2*B)) #    Equation du temps
    Lst = 15. * GMT # GMT-1  (pour la france)
    Lloc = longitude  # longitude du lieu
    
#Calcul de l'heure solaire    
    ts = heure_23 + (4 * (Lst - Lloc) + E) / 60.
    
#Angle horaire (Temps solaire en angle)
    omega = (ts - 12) * 15 * deg2rad
    #Position apparente du soleil
    
###Angle du zénith (theta_z)###
    phi = latitude * deg2rad # latitude
    thetaz = np.cos(phi) * np.cos(declinaison) * np.cos(omega) + np.sin(phi) * np.sin(declinaison)
    theta_z2 = np.arccos(thetaz)
    theta_z2[theta_z2>pi/2]=pi/2 #lorsque theta_z>90°, le soleil est couché
    theta_z =  theta_z2 * rad2deg
    
### Azimute solaire (gamma_s)###     
    gamma = (np.cos(theta_z * deg2rad) * np.sin(phi) - np.sin(declinaison)) / (np.sin(theta_z 
    * deg2rad) *np.cos(phi))
    gamma_s = np.sign(omega) * np.arccos(gamma) * rad2deg
    
### Angle d'incidence pour une surface orienté au sud (theta3_s)###
    gamma0 = 0 #azimute de la surface
    theta_s = np.cos(theta_z * deg2rad) * np.cos(beta*deg2rad) + np.sin(theta_z * 
    deg2rad) * np.sin(beta * deg2rad) * np.cos((gamma_s-gamma0)*deg2rad) #a verifier
    theta_s[theta_s<0]=0    
            
### Angle d'incidence pour une surface orienté au nord (theta3_n)###
    gamma0 = 180 #azimute de la surface
    theta_n = np.cos(theta_z * deg2rad) * np.cos(beta*deg2rad) + np.sin(theta_z * 
    deg2rad) * np.sin(beta * deg2rad) * np.cos((gamma_s-gamma0)*deg2rad) 
    theta_n[theta_n<0]=0
            
### Angle d'incidence pour une surface orienté à l'est (theta3_e)###
    gamma0 = -90 #azimute de la surface
    theta_e = np.cos(theta_z * deg2rad) * np.cos(beta*deg2rad) + np.sin(theta_z * 
    deg2rad) * np.sin(beta * deg2rad) * np.cos((gamma_s-gamma0)*deg2rad)
    theta_e[theta_e<0]=0

### Angle d'incidence pour une surface orienté à l'ouest (theta3_o)###
    gamma0 = 90 #azimute de la surface
    theta_o = np.cos(theta_z * deg2rad) * np.cos(beta*deg2rad) + np.sin(theta_z * 
    deg2rad) * np.sin(beta * deg2rad) * np.cos((gamma_s-gamma0)*deg2rad)   
    theta_o[theta_o<0]=0

    
    I_dir_s = I_dir * theta_s
    I_dir_n = I_dir * theta_n
    I_dir_e = I_dir * theta_e
    I_dir_o = I_dir * theta_o
    I_ref = rho * (I_dir * np.cos(theta_z)+I_dif) * (1 + np.cos(beta*deg2rad)) / 2.#Isotropique
    I_dif_surf = I_dif * (1 + np.cos(beta * deg2rad)) / 2. #isotropique
    I_ref[I_ref<0] = 0
    
    #I_moy = (I_dir_s+I_dir_n+I_dir_e+I_dir_o) / 4. + I_dif_surf + I_ref
    
    I_moy = np.zeros([8760*3600/delta,NB_dwe])
    for i in range(0,NB_dwe):
        temp = city['R_window_s'].values[i]*I_dir_s + city['R_window_n'].values[i]*I_dir_n + city['R_window_e'].values[i]*I_dir_e +\
        city['R_window_w'].values[i]*I_dir_o + I_dif_surf + I_ref
        I_moy[:,i] = tr.transformer(temp,3600,delta)
    
    I_dir_n = tr.transformer(I_dir_n,3600,delta)
    I_dir_s = tr.transformer(I_dir_s,3600,delta)
    I_dir_e = tr.transformer(I_dir_e,3600,delta)
    I_dir_o = tr.transformer(I_dir_o,3600,delta)
    I_dif_surf = tr.transformer(I_dif_surf,3600,delta)
    I_ref = tr.transformer(I_ref,3600,delta)
    I_dif = tr.transformer(I_dif,3600,delta)
    I_dir = tr.transformer(I_dir,3600,delta)
    
    return I_dir_s, I_dir_n, I_dir_e, I_dir_o, I_ref, I_dif_surf, I_moy, I_dir, I_dif
       