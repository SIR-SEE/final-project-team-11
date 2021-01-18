from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# Let’s run the basic SIR model

# describe the model
def deriv(y, t, N, beta, k, delta, m ):
    S, E, I, R, D = y
    dSdt = -beta  * S * I / N
    dEdt = (beta  * S * I / N - k * E) 
    dIdt = delta * E - k * I 
    dRdt = k * I 
    dDdt = m * I
    return  dSdt, dEdt, dIdt, dRdt, dDdt

# describe the parameters
N =  1000000 # population
m = 0.02    #percentage of the contaminated that dies    
delta = 0.2
P = 0.9 #the effect of restrictions in percent, 1 = no effect 0 = full effect
beta = P * 2.5 
k=1/6                   
S0, E0, I0, R0, D0= N-1, 1, 0, 0, 0 # initial conditions: one infected, rest susceptible


t = np.linspace(0, 99, 100) # Grid of time points (in days)
y0 = S0, E0, I0, R0, D0 # Initial conditions vector

    # Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, k, delta, m, ))
S, E, I, R , D= ret.T




def plotsir(t, S, E, I, R, D):
    f, ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
    ax.plot(t, E, 'y', alpha=0.7, linewidth=2, label='Exposed')
    ax.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
    ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')
    ax.plot(t, D , 'k', alpha=0.7, linewidth=2, label='Deceased')
    ax.set_xlabel('Time (days)')
    
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    #plt.savefig('C:/Users/albin/OneDrive/Dokument/GitHub/SIR_model_SEE070/SIR_model_SEE070/Lecture2dec.png')
    plt.show()
plotsir(t, S, E, I, R, D)
