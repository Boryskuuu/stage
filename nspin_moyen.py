# J'essaie de créer et applique l'opérateur évolution !

from qutip import * #analysis:ignore
import numpy as np
import matplotlib.pyplot as plt

### Fonctions à réutiliser

def idqbit(N):      #Crée matrice identité compatible format tenseurs
    l = []
    for i in range (N):
        l.append(identity(2))
    return(tensor(l))

N=50 # Nombre de spin auxiliaires
a=1
g=1
b=0.01

### Creation de H 


H=tensor(a*sigmax(),identity(2))    #Ha
H+=tensor(identity(2),N*b*sigmaz()) #Hs
H+=tensor(g*sigmax(),N*sigmax())    #Hint

## Création du vecteur d'état

psi=tensor(basis(2,1),basis(2,0))

### Évolution dans le temps

dt=0.0005
Nt = 50000

## Opérateur évolution de dt 

U=(-1j*H*dt).expm()    # fonction de qutip super pratique !

lp = [psi.ptrace(0)[0][0][0]]    
lm = [psi.ptrace(0)[1][0][1]]
lpa = [psi.ptrace(1)[0][0][0]]    
lma = [psi.ptrace(1)[1][0][1]]     

for i in range (Nt) :
    psi = U*psi
    lp.append(psi.ptrace(0)[0][0][0])   
    lm.append(psi.ptrace(0)[1][0][1])
    lpa.append(psi.ptrace(1)[0][0][0])
    lma.append(psi.ptrace(1)[1][0][1])

lt = np.linspace(0,(Nt+1)*dt,Nt+1)
plt.plot(lt,lp,label = "plus")
plt.plot(lt,lm,label = "moins")
plt.plot(lt,lpa,label = "plus aux moy")
plt.plot(lt,lma,label = "moins aux moy")
plt.grid()
plt.legend()


