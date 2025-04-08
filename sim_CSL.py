import numpy as np
import matplotlib.pyplot as plt

def gauss(sig, mu, a, x):
    g = a * np.exp(-(x - mu) ** 2 / (2 * sig ** 2))
    return g

def density(rc,phi,x,y) : 
    gaussians = gauss(rc, x[:,np.newaxis], np.abs(phi) ** 2 / (rc*np.sqrt(2 * np.pi)), y)
    return np.sum(gaussians,axis=0)

def Hcentre(N):     # Hamiltonien qui amène vers le centre
    H = np.zeros((N,N))
    for i in range (round(N/2)):
        H[i][i+1] = 1
    for i in range(round(N/2),N):
        H[i][i-1] = 1
    return(H)

def Hdrift(N):  # Hamiltonien qui amène vers la droite
    H = np.zeros((N,N))
    for i in range (N-1) :
        H[i][i+1] = 1
    H[N-1][N-1] = 1
    return(H)
        

def timeevol(phi,rc,gamma,x,dt,dx,N,H):
    schro = -1j*dt*np.matmul(H,phi)
    dens = density(rc,phi,x,x)
    dW = np.sqrt(dt) * np.random.normal(0, 1, size=N)
    wiener = dens * dW * phi 
    deterministic = - 0.5 * gamma * dens**2 * dt * phi
    # update + renormalize
    phi += wiener + deterministic + schro
    phi += schro
    phi /= np.sqrt(np.sum(np.abs(phi)**2) * dx)
    return(phi)
            
def doublegauss(sigma,center,a,x,dx) :
    phi = gauss(sigma, center, a, x) + gauss(sigma, -center, a, x)
    phi /= np.sqrt(np.sum(np.abs(phi)**2) * dx)  # Normalize
    return(phi)

def runvariance(sigma,center,a,x,dx,Nt,N,dt):
    phi = doublegauss(sigma,center,a,x,dx) 
    indvar = 0
    lvar = np.zeros(Nt)
    conv = 0
    stop = 0
    for i in range(Nt) :
        phi = timeevol(phi,rc,gamma,x,dt,dx,N)
        lvar[i] = np.var(x*np.abs(phi)**2)
        if np.abs(lvar[i]-lvar[i-1])<1e-8 :
            conv += 1
        else : 
            conv = 0
        if conv == 3 and stop == 0 :
            indvar[c] = (i-3)
            stop = 1
    return(indvar)

    
    
def runsingle(sigma,center,a,x,dx,Nt,N,dt):
    phi = doublegauss(sigma,center,a,x,dx).astype(complex)
    H = np.zeros((N,N))
    for i in range(Nt) :
        phi=timeevol(phi,rc,gamma,x,dt,dx,N,H)
        if i%50 == 0 :
            plt.plot(x, np.abs(phi) ** 2, label=f"t={i*dt:.3f}")            
    # Plot densité de proba 
    #plt.plot(x, np.abs(phi) ** 2, label="final")
    plt.xlabel("x")
    plt.ylabel("Probability Density")
    plt.title("|φ(x)|²")
    plt.legend()
    plt.grid()
    plt.show()



# Parameters
sigma = 1     # Initial width
N = 1000      # Number of spatial points
D = 25        # Max spatial range
a = 1.0  # Amplitude before normalization
rc = 1000
dx = 2 * D / N
x = np.linspace(-D, D, N)
lambda_csl = 1e-6
gamma = lambda_csl * (4 * np.pi * rc**2 ) ** (3 / 2)  # From CSL model
Nt = 500     # Number of time steps
dt = 0.01


# code pour run une seule simu (avec plots...)
center=6
runsingle(sigma,center,a,x,dx,Nt,N,dt)


# code pour run la variance (avec plots...)
# cmax = 60 
# cfreq = 60
# indvar = np.zeros(cmax)
# for c in range(cmax) :
#     center = c/cfreq    # Gaussian centers
#     nruns = 10
#     for run in range(nruns) : 
#         indvar[c] += runvariance(sigma,center,a,x,dx,Nt,N,dt)
#     indvar[c] /= nruns
# # Plot moment variance
# plt.plot(np.linspace(0,cmax/20,cmax),indvar)
# plt.xlabel("distance centre gaussiennes")
# plt.ylabel("étape du collapse")  
