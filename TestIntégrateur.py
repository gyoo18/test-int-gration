import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from math import *
import sys
from PolynômesPO import Poly

dt = 0.1            # delta temps en secondes
n_bonds = 20      # nombre de bonds à effectuer

k = 10      # Constante de rappel du ressort en N/m -> kg/s²
v0 = 0      # Vitesse initiale de la masse en m/s
x0 = 10     # Position initiale de la masse en m
m = 10      # Masse en kg de la masse

def calculer_sinusoïdale() -> list:

    global dt
    global n_bonds

    global k
    global v0
    global x0
    global m

    positions = [0]*n_bonds

    
    for i in range(n_bonds):
        positions[i] = x0*cos(sqrt(k/m)*i*dt)+((v0/sqrt(k/m)) * sin(sqrt(k/m)*i*dt))
        # positions[i] = 1/((i*dt)-2.0001)
    
    return positions

def intégration_euler() -> list:
    global dt
    global n_bonds

    global k
    global v0
    global x0
    global m

    positions = [0]*n_bonds
    positions[0] = x0

    x = x0
    v = v0

    for i in range(n_bonds-1):
        a = 2/(x*x*x)
        v += a*dt
        x += v*dt
        positions[i+1] = x

    return positions

def intégration_verletA() -> list:
    global dt
    global n_bonds

    global k
    global v0
    global x0
    global m

    positions = [0]*n_bonds
    positions[0] = x0

    x = x0
    xp = x0-(v0*dt)
    a = 0

    for i in range(n_bonds-1):
        a = -k*x/m
        xt = x
        x = 2*x - xp + a*dt*dt
        positions[i+1] = x
        xp = xt

    return positions

def intégration_verletB() -> list:
    global dt
    global n_bonds

    global k
    global v0
    global x0
    global m

    positions = [0]*n_bonds
    positions[0] = x0

    x = x0
    v = v0
    a = 0

    for i in range(n_bonds-1):
        v += 0.5*(-k*x/m)*(dt)
        x = x + v*dt
        a = -k*x/m
        v += 0.5*a*(dt)
        positions[i+1] = x

    return positions

def intégration_polynomiale4(n : int, s : int) -> list:
    global dt
    global n_bonds

    global k
    global v0
    global x0
    global m

    positions = [0]*n_bonds*s
    accél = [0]*n_bonds*s

    # Cette intégration a besoin des 3 derniers points, 
    # donc elle utilise Verlet-vitesse pour ces 3.

    x = x0
    xp = x0
    v = v0
    a = -k*x0/m

    positions[0] = x0
    accél[0] = a

    for i in range(n):
        for j in range(s):
            v += 0.5*(-k*x/m)*(dt/s)
            xp = x
            x = x + v*(dt/s)
            a = -k*x/m
            v += 0.5*a*(dt/s)
            positions[i*s+j+1] = x
            accél[i*s+j+1] = a
    
    # Intégration polynomiale

    for i in range(n_bonds-n-1):
        # coefficients = Poly([j for j in range(n+1)],accél[i*s+n*s-n:i*s+n*s+1])
        # print(coefficients)
        for k in range(s):
            coefficients = Poly([j for j in range(n+1)],accél[i*s+k+n*s-n:i*s+k+n*s+1])
            xt = x
            x = 2*x-xp + a*(dt/s)**2
            diviseur = 6
            for j in range(n//2):
                x += (1/diviseur)*coefficients[2*(j+1)]*(dt/s)**(2*j+4)
                diviseur *= (2*j+5)*(2*j+6)
            positions[i*s+k+n*s+1] = x
            xp = xt
            a = -k*x/m
            accél[i*s+k+n*s+1] = a
    
    posFin = [0]*n_bonds
    for i in range(n_bonds):
        posFin[i] = positions[i*s]
    
    return positions

def intégration_polynomiale_n(n : int):

    if type(n) != int or n%2 != 0 or n < 2:
        raise ValueError("[Intégration_polynomiale_n] le paramètre n doit être un entier pair et plus grand ou égal à 2.")

    global dt
    global n_bonds

    global k
    global v0
    global x0
    global m

    positions = [0]*n_bonds
    polynôme : list[list] = []
    for i in range(n//2-1):
        polynôme.append([0]*n_bonds)

    x = x0
    xp = x0
    v = v0
    a = -k*x0/m

    positions[0] = x0

    # xxxxxxx
    # vvvvvvv
    # aaaaaaa
    #  cccccc
    #4  bbbbb
    #    dddd
    #6    eee
    #      ff
    #8      g
    for i in range(n-1):
        v += 0.5*(-k*x/m)*(dt)
        xp = x
        x = x + v*dt
        a = -k*x/m
        v += 0.5*a*(dt)
        positions[i+1] = x
        if len(polynôme) != 0:
            polynôme[0][i+1] = a
        if i >= 2:
            for j in range(len(polynôme)-1):
                polynôme[j+1][i+1] = 2*( ( polynôme[i+1-3] - polynôme[i+1-2] )*( (i+1-2)*dt - (i+1-1)*dt ) - ( (i+1-3)*dt - (i+1-2)*dt )*( polynôme[i+1-2] - polynôme[i+1-1] ))/( ((i+1-3)*dt-(i+1-1)*dt)*((i+1-3)*dt-(i+1-2)*dt)*((i+1-2)*dt-(i+1-1)*dt) )
    
    # Intégration polynomiale

    for i in range(n_bonds-n):
        c = 0
        for j in range(len(polynôme)):
            polynôme[j][i+n] = 2*( -( polynôme[j][i+n-3] - polynôme[j][i+n-2] )*dt + ( polynôme[j][i+n-2] - polynôme[j][i+n-1] )*dt )/(-2*dt*dt*dt)
            c += ( 2/factorial( 2*(j+1) ) )*polynôme[j][i+(n-1)]*pow(dt,2*(j+1))
        a = -k*x/m
        xt = x
        x = 2*x - xp + a*dt*dt + c
        positions[i+n] = x
        xp = xt
    
    return positions



X = [0]*n_bonds
for i in range(len(X)):
    X[i] = i
    
vrai_fonction = calculer_sinusoïdale()
# euler = intégration_euler()
# verletA = intégration_verletA()
# verletB = intégration_verletB()
polynomiale4 = intégration_polynomiale4(2,1)
# polynomiale_n = intégration_polynomiale_n(4)

plt.plot(X,vrai_fonction, 'red')
# plt.plot(X,euler, 'green')
# plt.plot(X,verletA, 'blue')
# plt.plot(X,verletB, 'orange')
plt.plot([i/1 for i in range(n_bonds*1)],polynomiale4, 'cyan')
# plt.plot(X,polynomiale_n, 'purple')
# plt.gca().set_xlim([n_bonds-100, n_bonds])
# plt.gca().set_xlim([0, 10])
# plt.gca().set_ylim([-2, 2])
plt.show()