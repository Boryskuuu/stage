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
    li = [g*sigmax()]
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

U=(-1j*H*dt).expm()   #commande de base de Qutip !

lp = [psi.ptrace(0)[0][0][0]]    
lm = [psi.ptrace(0)[1][0][1]]
lpa = [psi.ptrace(5)[0][0][0]]    
lma = [psi.ptrace(5)[1][0][1]]     

for i in range (Nt) :
    psi = U*psi
    lp.append(psi.ptrace(0)[0][0][0])   
    lm.append(psi.ptrace(0)[1][0][1])
    lpa.append(psi.ptrace(5)[0][0][0])
    lma.append(psi.ptrace(5)[1][0][1])

lt = np.linspace(0,(Nt+1)*dt,Nt+1)
plt.plot(lt,lp,label = "plus")
plt.plot(lt,lm,label = "moins")
plt.plot(lt,lpa,label = "plus aux")
plt.plot(lt,lma,label = "moins aux")
plt.grid()
plt.legend()


