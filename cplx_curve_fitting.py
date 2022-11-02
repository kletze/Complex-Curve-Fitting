#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 17:02:19 2022
Complex curve fitting algorithm by E.C.Levy in IRE Transactions on Automatic Control AC-4 (1959)
@author: Stefan Kletzenbauer
@mail:   stefan.kletzenbauer@kern-microtechnik.com
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
# Funktion
# Input:
    # AMplitude, Winkel in Grad, zaehlerordnung, nennerordnung
def cplx_curve_fit( mag, phase, omega, ord_num, ord_den):
    realval = mag * np.cos(phase*np.pi/180)
    imagval = mag * np.sin(phase*np.pi/180)
    
# def cplx_curve_fit( cplx, omega, ord_num, ord_den):

#     realval = np.real(cplx)
#     imagval = np.imag(cplx)

    hmax = 2*ord_num+1

    h = 0

    lam = np.zeros(hmax)
    S = np.zeros(hmax)
    T = np.zeros(hmax)
    U = np.zeros(hmax)
    for h in range(0,hmax):

        lam[h] = np.sum(omega**h)
        S[h] = np.sum(omega**h * realval)
        T[h] = np.sum(omega**h * imagval)
        U[h] = np.sum(omega**h * (imagval**2 + realval**2))

    # lambda matrix

    ml = np.zeros((ord_num+1, ord_num+1))

    for i in range(ord_num+1):
        k = 1
        for j in range(ord_num+1):
            if not(j%2) and (j+(i%2) < ord_num+1 ):
                ml[i, j + (i%2)] = k * lam[i+j+(i%2)]
                k = -1*k

    print(lam)
    print("\n")
    print(ml)


    mden = np.zeros((ord_den+1, ord_den))
    # Matrix rechts hat n=ord_den Spalten und ord_den+1 Zeilen
    vzz = 1
    for z in range(ord_den+1): #Zeilen
        vzs = 1
        if z%2: #bei ungerader Zeilennummer:
            vzs = -1
            for s in range(ord_den):
                if s%2:
                    mden[z,s] = vzs * T[s+z+1]
                else:
                    mden[z,s] = vzs * S[s+z+1]
                    vzs = -1*vzs

        else:
            vzs = 1
            for s in range(ord_den): #Spalten bei gerader Zeilennummer
                if s%2: #bei ungerader Spaltennummer
                    mden[z,s] = vzs * S[s+z+1]
                    vzs = -1*vzs
                else:   #bei gerader Spaltennummer
                    mden[z,s] = vzs *T[s+z+1]


     #print("Mtrx: ")
     #print(mden)

     #Matrix unten
    munten = np.zeros((ord_den, ord_den+1))
    # Matrix rechts hat n=ord_den+1 Spalten und ord_den Zeilen
    vzz = 1
    for z in range(ord_den): #Zeilen
        vzs = 1
        if z%2: #bei ungerader Zeilennummer:
            vzs = 1
            for s in range(ord_den+1):
                if s%2:
                    munten[z,s] = vzs * T[s+z+1]
                    vzs = -1*vzs
                else:
                    munten[z,s] = vzs * S[s+z+1]


        else:
            vzs = 1
            for s in range(ord_den+1): #Spalten bei gerader Zeilennummer
                if s%2: #bei ungerader Spaltennummer
                    munten[z,s] = vzs * S[s+z+1]

                else:   #bei gerader Spaltennummer
                    munten[z,s] = vzs *T[s+z+1]
                    vzs = -1*vzs


    # U-Matrix:
    U_ger = U[0:len(U):2]
    mu = np.eye(len(U_ger))*U_ger
    muerg = mu[1:ord_num+1, 1:ord_num+1]

    LinkeM = np.concatenate((ml, munten))
    RechteM = np.concatenate((mden, muerg))

    sizeL = np.shape(LinkeM)
    sizeR = np.shape(RechteM)

    Mges = np.zeros((sizeL[0], sizeR[1]+sizeL[1]))
    for i in range(np.shape(LinkeM)[0]):
        Mges[i] = np.append(LinkeM[i], RechteM[i])

    C = np.zeros(ord_num+ord_den+1)
    for i in range(len(C)):
        if i <= ord_num:
            if i%2: #ungerade
                C[i] = T[i]
            else: # gerade
                C[i] = S[i]
        else:
            if (i-ord_num+1)%2:
                C[i] = U[i-ord_num]
            else:
                C[i]= 0

    print("Mges")
    print(Mges)
    N = np.linalg.solve(Mges, C)
    Numerator = N[0:ord_num+1]
    Denominator =  N[ord_num+1:ord_num+ord_den+1]
    Denominator = np.insert(Denominator, 0, 1.0)

    print(Numerator)
    print(Denominator)
    return Numerator, Denominator
