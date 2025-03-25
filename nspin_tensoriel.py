# J'essaie de créer et applique l'opérateur évolution !

from qutip import *
import numpy as np

N=10 # Nombre de spin auxiliaires
a=1
g=1
b=0.01

### Creation de H 

la=[a*sigmax()]     # Création de Ha
for i in range (N) :
    la.append(identity(2))
H=tensor(la)

for i in range (N) :
    ls = [identity(2)]
    li = [(g/N)*sigmax()]
    for j in range (N):
        if i==j :
            ls.append(b*sigmaz())
            li.append(sigmax())
        else : 
            ls.append(identity(2))
            li.append(identity(2))
    H+=tensor(li)+tensor(ls)

## Création du vecteur d'état

lp = [basis(2,1)]       # etat initial de l'atome
for i in range (N) :
    lp.append(basis(2,0))
psi=tensor(lp)

### Évolution dans le temps

dt=0.001
Nt = 50000

## Opérateur évolution de dt 

U = identity(psi.shape[0]) - 1j*dt*H

