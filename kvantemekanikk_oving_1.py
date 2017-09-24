import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

N = 500 #antall punkter 
hbar = 1
L = 20
res = 0.01
dx = L/(N-1)
Vmax = 0
xs = 5
sigma = 1
k = 20
m = 1
dt = res*(hbar/(hbar**2/(2*m*(dx**2))+Vmax))
E = hbar**2*k**2/(2*m)
w = E/hbar
vg = hbar*k/m
tmax = L/(2*vg)
print("Antall tiddsteg for å komme halveis i stedet: ", int(tmax/dt))
C = 1/(np.pi*sigma**2)**(1/4) #normeringskonstant
x = np.linspace(0,L,N) #first array
psiR = C*np.exp(-(x-xs)**2/(2*sigma**2))*np.cos(k*x) #startvalues for (x,0)
psiI = C*np.exp(-(x-xs)**2/(2*sigma**2))*np.sin(k*x-w*dt/2)
psiR[0] = 0
psiR[N-1] = 0
psiI[0] = 0
psiR[N-1] = 0
V = np.zeros(N)
V[int(N/2):int(N/2)+100]=np.linspace(E/2,1.2*E,100)
#her sender vi inn funksjonene, og returneres igjen etter ett tidssteg dvs dt
def step(V,dx,dt,psiR,psiI):
    psiI[1:N-1] = psiI[1:N-1]-dt*(
    (V[1:N-1]/hbar)*psiR[1:N-1]-(hbar/(2*m))*(psiR[2:N]-2*psiR[1:N-1]+psiR[0:N-2])/(dx)**2)
    #for å finne psiR må vi først oppdatere psiI!
    psiR[1:N-1] = psiR[1:N-1]+dt*((V[1:N-1]/hbar)*psiI[1:N-1]-
    (hbar/(2*m))*(psiI[2:N]-2*psiI[1:N-1]+psiI[0:N-2])/(dx)**2)
    return psiR, psiI
    
def normalize(dx,psiR,psiI):
    Psi = psiR**2+psiI**2
    return np.sum(Psi*dx)

for i in range(0,int(tmax/dt)):
    psiR,psiI = step(V,dx,dt,psiR,psiI)

a = normalize(dx,psiR,psiI)
print("Normering er: ", a)
plt.plot(x,psiR, linewidth = 0.5)
plt.plot(x,psiI, linewidth = 0.5)
plt.plot(x,psiR**2+psiI**2)
plt.plot(x,V/E, linewidth = 0.9)
plt.grid()
plt.show()

def f():
    return "Halla"
print(f)
"""
##animering

prob = [psiR**2+psiI**2]
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
line, = ax.plot([], [], lw=1)

for i in range(0,int(30*tmax/dt)):
    psiR,psiI = step(V,dx,dt,psiR,psiI)
    prob.append(psiR**2+psiI**2)
    
prob = prob[0:-1:100]

# Initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# Animation function which updates figure data.  This is called sequentially
def animate(i):
    line.set_data(x/L, prob[i])
    return line,

# Call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(prob), interval=1, blit=True)
plt.plot(x/L,V/E)
plt.grid()
plt.show()
"""