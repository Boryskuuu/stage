# CSL pour 1 spin
# J'essaie de reproduire les graphiques de l'article sur dCSL
# cf. La figure 1 
# etendre à la dCSL après ? / plusieurs particules ?

import numpy as np
import matplotlib.pyplot as plt

def gauss(sig,mu,a,x) :
    g = a * np.exp(-(x-mu)**2/(2*sig**2))
    return(g)

sigma = 1 # variance initiale...
center = 4 # Centre des gausiennes

# Discrétisation de l'espace 

N = 1000 # nombres de points espace
D = 25 # distance max à 0  

x = np.linspace(-D,D,N)

# génération de l'état initial 

a=1/np.sqrt(2*sigma*np.sqrt(2*np.pi))   # car normalise pour norme

phi= gauss(sigma,center,a,x) + gauss(sigma,-center,a,x)

# Evolution dans le temps 
# temps = t/lambda

Nt = 1
dt=1/Nt
temps = np.linspace(0,1,Nt)

# On exprime en lambda=1 et rc=1 

gamma = (4*np.pi)**(3/2)    # voir formule papier 
m = 1

def smeared_mass(phi,x,m): # a complexifier pour N parts...
    M = 0   ## passer en VECTORIEL !!!!
    for k in range(len(phi)) :
        M += gauss(1,x,m * phi[k]**2 / (2*np.pi)**(3/2),k/len(phi))
    return(M)

plt.plot(x,phi**2,"g")  # phi**2 pour norme

moyenne = np.mean(phi)

## le processus de Wiener devrait s'appliquer de manière homogène sur tout le vecteur....


for i in range(Nt):
    dW = np.random.normal(0, 0.1, size=N)
    for k in range(N):
        phi[k] += - np.sqrt(gamma) / m * moyenne  * dW[k]
        phi[k] += - dt**2 * gamma / (2 * m**2) * moyenne**2

# for i in range(Nt):
#     dW = np.random.normal(0, 1, size=N)
#     save=smeared_mass(phi,x,m)
#     phi += sum(smeared_mass(phi,x,m)-m)
#     phi += sum(smeared_mass(phi,x,m)-m)
#     # phi += np.sqrt(gamma) / m * sum(smeared_mass(phi,x,m)-m)  * dW
#     # phi += - dt**2 * gamma / (2 * m**2) * sum(smeared_mass(phi,x,m)-m)

plt.plot(x,phi**2,"r")  # phi**2 pour norme
plt.grid()
plt.show()

