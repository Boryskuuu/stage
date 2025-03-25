# J'essaie de créer et applique l'opérateur évolution !

from qutip import *
import numpy as np

N=10
U=identity(2)
for i in range (N) :
    U=tensor(U,identity(2))
print(U)
