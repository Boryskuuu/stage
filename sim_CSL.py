import numpy as np
import matplotlib.pyplot as plt


def gauss(sig, mu, a, x):
    g = a * np.exp(-(x - mu) ** 2 / (2 * sig ** 2))
    return g

def density(rc,phi,x,y) : 
    gaussians = gauss(rc, x[:,np.newaxis], np.abs(phi) ** 2 / (rc*np.sqrt(2 * np.pi)), y)
    return np.sum(gaussians,axis=0)


def smeared_mass(phi, x, y, m, rc):
    gaussians = gauss(rc, x[:,np.newaxis], m * np.abs(phi) ** 2 / (rc*np.sqrt(2 * np.pi)), y)
    return sum(gaussians)


# Parameters
sigma = 1     # Initial width
N = 1000      # Number of spatial points
D = 25        # Max spatial range
rc = 1000
dx = 2 * D / N
x = np.linspace(-D, D, N)
lambda_csl = 0.01
H = 0 
cmax = 60 
indvar = np.zeros(cmax)

for c in range(cmax) :
    # Generate initial wave function
    center = c/20    # Gaussian centers
    a = 1.0  # Amplitude before normalization
    nruns = 10
    for run in range(nruns) : 
        phi = gauss(sigma, center, a, x) + gauss(sigma, -center, a, x)
        phi /= np.sqrt(np.sum(np.abs(phi)**2) * dx)  # Normalize
        
        # Time evolution
        Nt = 100     # Number of time steps
        dt = 1 / Nt
        gamma = lambda_csl * (4 * np.pi * rc**2 ) ** (3 / 2)  # From CSL model
        m = 1
        
        # plt.plot(x, np.abs(phi) ** 2, label="t=0")
        
        lvar = np.zeros(Nt)
        conv = 0
        stop = 0
        
        for i in range(Nt) :
        #    schro = -1j*H*dt*phi
            dens = density(rc,phi,x,x)
            dW = np.sqrt(dt) * np.random.normal(0, 1, size=N)
            wiener = dens * dW * phi 
            deterministic = - 0.5 * gamma * dens**2 * dt * phi
            
            # update + renormalize
            phi += wiener + deterministic
            phi /= np.sqrt(np.sum(np.abs(phi)**2) * dx)
            # if i < 16 and i%2 == 0 :
            #     plt.plot(x, np.abs(phi) ** 2, label=f"t={i*dt:.3f}")
            lvar[i] = np.var(x*np.abs(phi)**2)
            if np.abs(lvar[i]-lvar[i-1])<1e-8 :
                conv += 1
            else : 
                conv = 0
            if conv == 3 and stop == 0 :
                indvar[c] += (i-3)/nruns
                stop = 1
        
# Plot moment variance
plt.plot(np.linspace(0,cmax/20,cmax),indvar)
plt.xlabel("distance centre gaussiennes")
plt.ylabel("étape du collapse")        
        

#plt.plot(lvar)

# # Plot densité de proba 
# plt.plot(x, np.abs(phi) ** 2, label="final")
# plt.xlabel("x")
# plt.ylabel("Probability Density")
# plt.title("|φ(x)|²")
# plt.legend()
# plt.grid()
# plt.show()
