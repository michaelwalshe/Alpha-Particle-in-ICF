# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 08:54:57 2019

@author: Michael
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 00:10:22 2019
Calculate the stopping power of a plasma using T-matrix model from Edie PhD thesis
@author: Michael
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt

#Gives you the coefficient dependent on degeneracy param x (based on gamma and charge) and which coefficient i you want
def coef(i, x):
    if i == 1:
        a = m.exp(1.78 * x - 1.01)
    elif i==2:
        a = m.exp(2.561 * x - 1.15)
    elif i==3:
        a = m.exp(3.141 * x - 2.41)
    elif i==4 and x > 0.548:
        a = m.exp(1.78 * x - 1.01)
    elif i==4 and (x<0.548 or x>-2.649):
        a = m.exp(1.78 * x - 1.01)
    elif i==4 and x < -2.649:
        a = m.exp(1.78 * x - 1.01)
    elif i==5 and x > -3.13:
        a = m.exp(1.78 * x - 1.01)
    elif i==5 and (x<0.548 or x<-3.13):
        a = m.exp(1.78 * x - 1.01)
    elif i==6:
        a = m.exp(1.78 * x - 1.01)
    elif i==7 and x>1.83:
        a = m.exp(1.78 * x - 1.01)
    elif i==7 and x<1.83:
        a = m.exp(1.78 * x - 1.01)
    elif i==8:
        a = m.exp(1.78 * x - 1.01)
    elif i==9:
        a = m.exp(1.78 * x - 1.01)
    return(a)

#Calculates dE/dx for energy of alpha particle E
def dE(E):
    T = 10 ** 7 #temp in K
    n_e = 10 ** 25 #density of electrons (&ions*2) in cm^-3
    m_alpha = 6.644E-24 #mass of alpha in kg
    m_e = 9.1E-28 #mass of e- in g
    Z_b = 2 #charge of alpha
    e = 4.8E-10 #charge of e-
    k = 1.38E-16 #boltzmann
    hbar = 1.055E-27 #planck
    a_B = 5.29E-9 #assumed bohr radius but not sure
    Gamma = e ** 2 * (m.pow(((4 * m.pi * n_e) / 3), 1/3) / (k * T)) #Gamma_ee coupling param from edie
    v_th = m.sqrt((k * T) / m_e) #thermal velocity of electrons
    omega = m.sqrt((4 * m.pi * n_e * e ** 2) / m_e) #electron plasma frequency
    c1 = (2 * m_e * v_th ** 2) / (hbar * omega)
    c2 = (3 * a_B * (k * T) ** 3) / (Z_b ** 2 * e ** 2 * hbar ** 2 * omega ** 2)
    vel = m.sqrt(2 * E / m_alpha) #velocity of alpha particles
    x = m.log(Z_b * m.pow(Gamma, 1.5)) #x_0 argument
    vbar = vel / v_th
    dEdxtop = (coef(1, x) * vbar + coef(2, x) * vbar ** 2 + coef(3, x) * vbar ** 4 * m.log(c1 * vbar ** 2 + 1))
    dEdxbot = 1 + coef(5, x) * vbar + coef(6, x) * vbar ** 2 + coef(7, x) * vbar ** 3 + coef(8, x) * vbar ** 4 + coef(9, x) * vbar ** 5 + coef(4, x) * c2 * vbar ** 7
    dEdx = (dEdxtop / dEdxbot) * ((3 * (k * T) ** 2) / (e ** 2)) #not sure I need the second bit
    return(dEdx)

if __name__ == "__main__":
    T_c = 10 ** 7
    E_MeV = np.linspace(0.01, 6.0, 100)
    E_erg = E_MeV * 10 ** 6 * 1.6E-19 * 10 ** -7
    T_mat_dEdx = np.zeros_like(E_MeV)
    for i in range(100):
        T_mat_dEdx[i] = dE(E_erg[i]) / 1.6E-8
    print(E_MeV)
    print(T_mat_dEdx)
    plt.plot(E_MeV, T_mat_dEdx, label='T-matrix')
    plt.xlabel('E[MeV]')
    plt.ylabel('-dE/dx [MeV / micrometer]')
    plt.savefig('t_matrix.png')
    plt.show()
