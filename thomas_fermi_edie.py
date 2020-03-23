# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 22:39:13 2019
Calculate the stopping power of a plasma using Thomas-Fermi model from Edie PhD thesis
Almost definitely shouldn't be seperate fcts but that's how I started that's how I'll finish
@author: Michael
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt

#Calculates the chemical potential for a certain temperature (using Ichimaru approx. from stat. physics of plasma)
def fermi_int(T):
    n_e = 10 ** 25 #dens of e- in cm^-3
    hbar = 1.054E-27 #planck constant in cgs
    m_e = 9.1E-31 #mass of e- in g
    E_F = ((hbar ** 2)/(2 * m_e)) * m.pow((3 * m.pi ** 2 * n_e), 2/3)
    theta = T / E_F
    I_half = (2/3) * m.pow(theta, -1.5)
    A = 0.25954
    B = 0.072
    b = 0.858
    const = (A * m.pow(theta, -(b + 1)) + B * m.pow(theta, -0.5*(b + 1))) / (1 + A * m.pow(theta, -b))
    mu = T * (1.5 * m.log(theta) + 4 *  m.log(4/(3 * m.sqrt(m.pi))) + const)
    return (mu, I_half)

#This calculates the arg. x_0 for our DT plasma with alpha particle beam
#needs energy of beam epsilon_b, temp of plasma
def x_0(eps_b, T_c):
    m_b = 6.644E-24 #mass of alpha particle in g
    m_c = 4.15E-24 #avg mass of D & T in g
    k = 1.38E-16 #boltzmann in cgs
    fermi = fermi_int(T_c)
    eta = fermi[0] / (k * T_c)
    I = fermi[1]
    a1 = m.sqrt((m_c * eps_b) / (m_b * k * T_c))
    a2 = m.sqrt(m.pi) / (2 * I * (1 + m.exp(-eta)))
    a3 = m.pow(a2, 1/3)
    x = a1 * a3
    print(fermi)
    return(x)

#Heaviside step fct for x (0 unless x positive)
def erf(x):
    if x < 0:
        return 0
    else:
        return 1

##Calculates Pade approximation for coulomb logarithm
def pade(eps_b, T_c):
    m_e = 9.1E-28
    k = 1.38E-16
    v_e = m.sqrt((k * T_c) / m_e)
    n_e = 10 ** 25
    e = 4.8E-10
    hbar = 1.055E-27
    omega = m.sqrt((4 * m.pi * n_e * e ** 2) / m_e) #electron plasma frequency
    x = x_0(eps_b, T_c)
    b1 = (2 * m_e * v_e ** 2) / (hbar * omega)
    b2 = (0.321 + 0.259 * x ** 2 + 0.0707 * x ** 4 + 0.05 * x ** 6) / (1 + 0.13 * x ** 2 + 0.05 * x **4)
    return(b1 * b2)

#calc dE/dx
def dE(eps_b, T_c):
    Z_b = 2
    e = 4.8E-10
    m_b = 6.644E-24
    m_c = 4.15E-24
    v_b = m.sqrt((2 * eps_b) / m_b)
    x = x_0(eps_b, T_c)
    c1 = (4 * m.pi * Z_b ** 2 * e ** 4) / (m_c * v_b ** 2)
    c2 = m.log(pade(eps_b, T_c)) * (erf(x) - (1 + m_c/m_b) * ((2/m.sqrt(m.pi)) * x * m.exp(-x ** 2)))
    print(c1)
    return(- c1 * c2)

if __name__ == "__main__":
    T_c = 10 ** 7 #temp in K
    E_MeV = np.linspace(0.01, 4.0, 50)
    E_erg = E_MeV * 10 ** 6 * 1.6E-19 * 10 ** -7
    dEdx = np.zeros_like(E_MeV)
    for i in range(50):
        dEdx[i] = -1 * dE(E_erg[i], T_c) / 1.6E-8
    print(dEdx)
    plt.plot(E_MeV, dEdx)
    plt.xlabel('E[MeV]')
    plt.ylabel('-dE/dx [MeV / micrometer]')
    plt.savefig('thomas_fermi.png')
    plt.show()
