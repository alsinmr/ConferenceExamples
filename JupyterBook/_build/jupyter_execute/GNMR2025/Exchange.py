#!/usr/bin/env python
# coding: utf-8

# # Chemical Exchange

# Chemical exchange is an important phenomenon in magnetic resonance, where a chemical or physical change in a molecule gives rise to a modulation in the chemical shift for a given spin. Depending on the rate of exchange between two or more chemical shifts, we may observe the effects of exchange using a number of different experiments.
# 
# The following tutorial will introduce this concept and investigate effects of various parameters on several different experiments. We'll use SLEEPY [http://sleepy-nmr.org](http://sleepy-nmr.org) to support our investigations, so take some time to experiment with the software as you go.
# 
# The basic formulas for calculating the influence of exchange on a 1D spectrum are the Bloch-McConnell equations:
# 
# H.M. McConnell. [*J. Chem. Phys.*](https://doi.org/10.1063/1.1744152), **1958**, 28, 430-431.

# In[1]:


# Make sure SLEEPY is installed and on the path
import sys
import os
if 'google.colab' in sys.modules: #Colab
    get_ipython().system('pip install sleepy-nmr')
elif 'USER' in os.environ and os.environ['USER']=='jovyan': #Binder
    get_ipython().system('pip install sleepy-nmr')
else:
    os.chdir('../../..')


# ## Setup

# In[2]:


import SLEEPY as sl
import numpy as np
import matplotlib.pyplot as plt
sl.Defaults['verbose']=False


# ## Simulation functions
# The hidden cell below (visibile in Colab) defines 5 different experiments that we'll apply in this tutorial:
# 
# * `OneD` : Simple, 1D spectrum
# * `EXSY` : Exchange Spectroscopy experiment. Set buildup=False to see a 2D spectrum. Set buildup=True to see time-dependent buildup of the peaks over multiple EXSY experiments
# * `CEST` : Chemical Exchange Saturation Transfer experiment
# * `R1p`  : $R_{1\rho}$ relaxation (relaxation under a spin-lock). 6 field strengths applied
# * `CPMG` : Carr-Purcell-Meiboom-Gill pulse train. 4 frequencies applied (matches the first four field strengths of $R_{1\rho}$ experiment

# In[101]:


def OneD(L,ax=None):
    seq=L.Sequence(Dt=1/4000)

    rho=sl.Rho(rho0='S0x',detect='S0p')
    rho.DetProp(seq,n=4096)
    ax=rho.plot(FT=True,ax=ax)
    return ax
    
def EXSY(L,dly=1,buildup=True,ax=None):
    dly0=np.linspace(0,dly,20) if buildup else [dly]
    if buildup:I00,I01,I10,I11=[],[],[],[]
    
    p1=1/(1+L.kex[1,0]/L.kex[0,1])
    
    
    seq=L.Sequence(Dt=1/2000)
    rho=sl.Rho(rho0='S0x',detect='S0p')
    v1=50000
    tpi2=1/v1/4
    
    for dly in dly0:
        t=[0,tpi2,dly+tpi2,dly+2*tpi2]
        seqX=L.Sequence().add_channel(0,t=t,v1=[v1,0,v1],phase=[-np.pi/2,0,np.pi/2])
        seqY=L.Sequence().add_channel(0,t=t,v1=[v1,0,v1],phase=[0,0,np.pi/2])
        twoD=sl.Tools.TwoD_Builder(rho=rho,seq_dir=seq,seq_in=seq,seq_trX=seqX,seq_trY=seqY)
        twoD(64,32).proc()
        if buildup:
            I11.append(twoD.Sreal.real[:64][:,:32].sum())
            I10.append(twoD.Sreal.real[:64][:,32:].sum())
            I01.append(twoD.Sreal.real[64:][:,:32].sum())
            I00.append(twoD.Sreal.real[64:][:,32:].sum())
    
    if buildup:
        I00,I01,I10,I11=[np.array(I)/I00[0]*p1 for I in [I00,I01,I10,I11]]
        if ax is None:ax=plt.subplots()[1]
        ax.plot(dly0*1e3,I00,label=r'$0\rightarrow 0$')
        ax.plot(dly0*1e3,I01,label=r'$0\rightarrow 1$')
        ax.plot(dly0*1e3,I10,label=r'$1\rightarrow 0$')
        ax.plot(dly0*1e3,I11,label=r'$1\rightarrow 1$')
        ax.set_xlabel('t / ms')
        ax.legend()
    else:
        if ax is None:
            ax=twoD.plot()
            ax.figure.set_size_inches([9,9])
        else:
            twoD.plot(ax=ax)    
        ax.set_zlim([0,twoD.Sreal.real.max()*1.2])
    return ax

def CEST(L,v1=25,dly=1,ax=None):
    rho=sl.Rho(rho0='S0z',detect='S0p')
    seq=L.Sequence(Dt=1/2000)
    t=[0,dly,dly+2.5e-6]
    sat=L.Sequence().add_channel(0,t=t,v1=[v1,1e5],phase=[0,np.pi/2])
    rho,seq,sat=rho.ReducedSetup(seq,sat)

    I=[]
    voff0=np.linspace(-2000,2000,200)
    for voff in voff0:
        rho.clear(data_only=True)
        t=[0,dly,dly+2.5e-6]
        sat.add_channel(0,t=t,v1=[v1,1e5],voff=[voff,0],phase=[0,np.pi/2])
        (sat*rho).DetProp(seq,n=128)
        I.append(rho.FT[0][125:225].real.sum())
        
    if ax is None:
        ax=plt.subplots(1,1)[1]
    ax.plot(voff0/1e3,I)
    ax.set_xlabel(r'$\nu$ / kHz')
    ax.set_xlim([voff0[-1]/1e3,voff0[0]/1e3])
    return ax

def R1p(L,dly=0.25,ax=None):
    v10=np.logspace(np.log10(2000),np.log10(40000),6)
    rho=sl.Rho(rho0='S0x',detect='S0p')
    leg=[]
    for v1 in v10:
        seq=L.Sequence(Dt=dly/500).add_channel(0,v1=v1)
        rho.clear()
        rho.DetProp(seq,n=500)
        ax=rho.plot(ax=ax)
        leg.append(rf'{v1/1e3:.1f} kHz')
    sc=L.expsys.Peq[0]/2 if L.Peq else 0.5
    ax.set_ylim([-0.05*sc,sc])
    ax.legend(leg)
    return ax

def CPMG(L,dly=0.25,ax=None):
    f0=np.logspace(np.log10(2000),np.log10(40000),6)
    v1=50000
    pi2=1/v1/4
    dly0=1/f0/4
    n0=(dly/dly0/8).astype(int)
    rho=sl.Rho(rho0='S0x',detect='S0p')
    leg=[]
    for f,dly,n in zip(f0,dly0,n0):
        t=[0,dly-pi2,dly+pi2,3*dly-pi2,3*dly+pi2,5*dly-pi2,5*dly+pi2,7*dly-pi2,7*dly+pi2,8*dly]
        seq=L.Sequence().add_channel(0,v1=[0,v1,0,v1,0,v1,0,v1,0],t=t,phase=[0,0,0,np.pi/2,0,0,0,np.pi/2,0])
        rho.clear()
        rho.DetProp(seq,n=n)
        ax=rho.plot(ax=ax)
        leg.append(rf'{f/1e3:.1f} kHz')
    sc=L.expsys.Peq[0]/2 if L.Peq else 0.5
    ax.set_ylim([-0.05*sc,sc])
    ax.legend(leg)
    return ax

pars=[(.01,.95),(1e-5,0.75),(.1,0.35),(1e-3,.6)]


# ## Chemical Exchange in 1D

# The most basic experiment for observing chemical exchange is a simple 1D experiment, where we exchange between chemical shifts. In this experiment, evolution of a coherence under the chemical shift depends on the state of the system, for example:
# 
# $$
# \begin{eqnarray}
# \frac{d}{dt}\langle\hat{S}^+\rangle_1^+&=&-i\Omega_1\langle\hat{S}^+\rangle_2 \:\mathrm{ : state 1} \\
# \frac{d}{dt}\langle\hat{S}^+\rangle_2^+&=&-i\Omega_2\langle\hat{S}^+\rangle_2 \:\mathrm{ : state 2}
# \end{eqnarray}
# $$

# In[4]:


ex0=sl.ExpSys(v0H=500,Nucs='1H')    #Experimental system at 500 MHz with one 1H
ex0.set_inter('CS',i=0,Hz=500)     #Add a chemical shift to the system
ex1=ex0.copy()                      #Copy system for second system in exchange 
ex1.set_inter('CS',i=0,Hz=-500)    #Add a chemical shift to the second system
L=sl.Liouvillian(ex0,ex1)               #Generate a Liouvillian for the system
seq=L.Sequence(Dt=1/(10000))       #Sequence with a 25 microsecond length

rho=sl.Rho(rho0='1Hx',detect='1Hp')
rho.DetProp(seq,n=100)
_=rho.plot()


# As we see, the system simply oscillates at the given input frequencies. Below, we can check what that frequency actually is, by taking the Fourier transform of the above signal, where we see that indeed it oscillates at +500 and -500 Hz.

# In[5]:


ax=rho.plot(FT=True)
_=ax.set_xlim([3,-3])


# However, we may couple the two states together with an exchange matrix, yielding broadening in the spectrum, or if fast enough, coalescence of the peaks.

# In[6]:


L.kex=sl.Tools.twoSite_kex(tc=5e-4)
rho=sl.Rho(rho0='1Hx',detect='1Hp')
rho.DetProp(seq,n=100)
_=rho.plot(FT=True)


# ### Exercise 1.1:
# Typing any of the above objects at the command line will produce a description of the object. Objects also have a plot function (`ex0`,`ex1` have a `plot_inter` function, although this is not so informative for isotropic interactions). Try this to get a feeling for what the critical settings are in this simulation. What is the role of each object in the simulation?

# In[7]:


# Use this cell to investigate the different objects


# ```{toggle}
# * ExpSys (ex0,ex1) : defines spins in the magnetic field, spins in the system, and interactions
# 
# * Liouvillian (L) : Contains each of the ExpSys objects and the Liouvillian matrix, as well as an exchange matrix (currently set to zeros).
# 
# * Sequence (seq) : Contains definitions of the field strengths, and defines length of propagation steps.
# 
# * Rho (rho) : Contains the initial state of the system and the detection matrix, as well as the signal.
# ```

# Now, assume that the system is in one of two states, with chemical shift frequencies $\Omega_1$ and $\Omega_2$. Depending on which state we're in, we have either
# 
# $$
# \begin{eqnarray}
# \frac{d}{dt}\langle\hat{S}^+\rangle&=&-i\Omega_1\langle\hat{S}^+\rangle \\
# \frac{d}{dt}\langle\hat{S}^+\rangle&=&-i\Omega_2\langle\hat{S}^+\rangle
# \end{eqnarray}
# $$
# 
# Then, if we hop randomly from state 1 to 2 and vice versa, with rate constants $k_{12}$ and $k_{21}$, we have the following exchange matrix.
# 
# $$
# \mathbf{k}=
# \begin{pmatrix}
# k_{12} & k_{21} \\ k_{12} &- k_{21}
# \end{pmatrix}
# $$
# 
# This can be combined into a single matrix including the oscillation, which yields the following differential equation:
# 
# $$
# \frac{d}{dt}
# \begin{pmatrix}
# \langle\hat{S}^+_1\rangle \\ \langle\hat{S}^+_1\rangle
# \end{pmatrix}
# =
# \begin{pmatrix}
# -i\Omega_1 - k_{12} & k_{21} \\ k_{12} & -i\Omega_2 - k_{21}
# \end{pmatrix}
# \cdot
# \begin{pmatrix}
# \langle\hat{S}^+_1\rangle \\ \langle\hat{S}^+_1\rangle
# \end{pmatrix}
# $$
# 
# This matrix indicates that if we are in state 1, the phase of $\hat{S}_1^+$ will oscillate with frequency of $\Omega_1$, but it will also diffuse with rate $k_{12}$ to state 2, where it will oscillate with frequency $\Omega_2$ and also potentially diffuse back to state one (with rate $k_{21}$).

# #### Some hints on 2x2 exchange matrices
# A 2x2 exchange matrix has two eigenvalues. One must always be zero, corresponding to the equilibrium of the exchange. Any evolution of the matrix occurs with a rate constant given by the other eigenvalue. This eigenvalue is equal to $-(k_{12}+k_{21})$. We therefore often define $k=(k_{12}+k_{21})/2$, and $\tau_c=1/(2k)$.

# ### Exercise 1.2:
# To understand exchange analytically, we often make some simplifying assumptions about the terms. First, we set $\Omega_1=\Delta\Omega/2$ and $\Omega_2=-\Delta\Omega/2$. We also assume symmetric exchange, such that $k_{12}=k_{21}=k$. This yields the following exchange matrix:
# 
# $$
# \mathbf{k}=
# \begin{pmatrix}
# -i\Delta\Omega/2 - k & k \\ k & i\Delta\Omega/2 - k
# \end{pmatrix}
# $$
# 
# Find the eigenvalues of this matrix, where each eigenvalue represents the evolution of a signal in the NMR spectrum. We won't solve for the eigenvectors here, but these also play a role in the contributions of each signal.
# 
# Under what conditions do we observe two separated chemical shifts? What is $T_2$ in this case? What is $T_2$ if there is only one chemical shift?
# 
# Hint: $T_2$ is the negative inverse of the real part of the eigenvalues, whereas the imaginary part yields oscillation frequencies
# 
# The answer is given in multiple steps, so if you're not sure how to start, go ahead and open the first answer.

# #### Part 1
# ```{toggle}
# Eigenvalues are found by calculating the determinant of the matrix minus $\lambda\cdot\mathbb{1}$, setting it to zero, and solving for the values of $\lambda$ ($\mathbb{1}$ is an identity matrix)
# 
# $$
# 0=\det
# \begin{pmatrix}
# -i\Delta\Omega/2-k-\lambda & k \\ k & i\Delta\Omega/2-k-\lambda
# \end{pmatrix}
# $$
# ```

# #### Part 2
# ```{toggle}
# The determinant for a 2x2 matrix is calculated by taking the (0,0) element times the (1,1) element minus the (0,1) element times the (1,0) element.
# 
# $$
# \begin{eqnarray}
# 0&=&(-i\Delta\Omega/2-k-\lambda)\cdot(i\Delta\Omega/2-k-\lambda)-k^2 \\
# 0&=&\lambda^2+2k\lambda+\frac{(\Delta\Omega)^2}{4}
# \end{eqnarray}
# $$
# ```

# #### Part 3
# ```{toggle}
# We can solve for the eigenvalues using the quadratic formula.
# 
# $$
# \begin{eqnarray}
# 0&=&ax^2+bx+c \\
# x&=&\frac{-b\pm\sqrt{b^2-4ac}}{2a}
# \end{eqnarray}
# $$
# 
# Then,
# 
# $$
# \lambda=\frac{-2k\pm\sqrt{4k^2-(\Delta\Omega)^2}}{2}
# $$
# ```

# #### Part 4
# ```{toggle}
# The real part of the eigenvalues is the negative inverse of $T_2$, whereas the imaginary part is the frequency. When $2k>\Delta\Omega$, then the eigenvalues are strictly real. This means that there is only a single frequency (at 0) in the spectrum, but a complex lineshape with two decay rates. On the other hand, if $2k<\Delta\Omega$, then two distinct frequencies emerge. The real part is then only $-k$, so $T_2=1/k$.
# 
# $$
# \begin{matrix}
# 2k>>\Delta\Omega & T_2=\infty,1/(2k) & \Omega=0 \\
# 2k>\Delta\Omega & T_2=1/(k\mp\sqrt{k^2-(\Delta\Omega/2)^2}) & \Omega=0 \\
# 2k=\Delta\Omega & T_2=1/k & \Omega=0 \\
# 2k<\Delta\Omega & T_2=1/k & \Omega=\pm\sqrt{(\Delta\Omega/2)^2-k^2} \\
# 2k<<\Delta\Omega & T_2=1/k & \Omega=\pm\Delta\Omega/2
# \end{matrix}
# $$
# Note that when $k>\Delta\Omega$, there are two peaks with the same apparent chemical shift. As $k$ gets large, the contribution of $T_2=1/(2k)$ becomes smaller, such that the peak becomes increasingly narrow, rather than wide.
# ```

# ### Exercise 1.3
# Above, we set up a simulation with two different chemical shifts in exchange ($\pm$500 Hz, so $\Delta\Omega=2*\pi*1000$ rad/s). Adjust $k$ to the different conditions in the solution above to observe the effect.

# In[8]:


k=0
L.kex=[[-k,k],[k,-k]]

_=OneD(L)


# In[10]:


fig=plt.figure(figsize=[12,6])
ax=[fig.add_subplot(2,3,k+1) for k in range(5)]
DeltaCS=L.H[0].expsys[0]['Hz']-L.H[1].expsys[0]['Hz']  #Change in chemical shift (Hz)
sqrt0=2*np.pi*DeltaCS/2   #This is where the square root is zero
k0=[50*sqrt0,2*sqrt0,sqrt0,sqrt0/2,sqrt0/50]
titles=[rf'$2k{sign}\Delta\Omega$' for sign in ['\gg','>','=','<','\ll']]

for a,k,title in zip(ax,k0,titles):
    L.kex=[[-k,k],[k,-k]]
    OneD(L,ax=a)
    a.set_title(title)
_=fig.tight_layout()


# The full Liouville matrix for the above system is an 8x8 matrix, but most terms do not actually enter into the calculation, reducing to the 2x2 matrix already discussed. SLEEPY will automatically pick out the required terms. We may then adjust the exchange matrix and extract the eigenvalues of the Liouville matrix for the conditions above, in order to find $T_2$ and the resonance frequencies.

# In[16]:


rho=sl.Rho('S0x','S0p')
rho_r,_=rho.ReducedSetup(L.U(.1))  
L_r=rho_r.L   # The reduced Liouvillian is not returned directly
for k in k0:
    L.kex=[[-k,k],[k,-k]]
    Lambda=np.linalg.eig(L_r[0].L(0))[0]
    print(f'k={k:.1} s^-1, '+\
          f'T2={-1e3/Lambda[0].real:5.2f} ms, T2={-1e3/Lambda[1].real:5.2f} ms,'+\
          f'f={-Lambda[0].imag/2/np.pi:6.1f} Hz, f={-Lambda[1].imag/2/np.pi:6.1f} Hz')


# Note that we have assumed symmetric exchange. However, similar effects occur for asymmetric exchange. These are simulated below for the same conditions. The `twoSite_kex` tool builds an exchange matrix for us, defined by a correlation time $1/(2k)$ and a population for the first state.

# In[17]:


fig=plt.figure(figsize=[12,6])
ax=[fig.add_subplot(2,3,k+1) for k in range(5)]
DeltaCS=L.H[0].expsys[0]['Hz']-L.H[1].expsys[0]['Hz']  #Change in chemical shift (Hz)
sqrt0=2*np.pi*DeltaCS/2   #This is where the square root is zero
k0=[50*sqrt0,2*sqrt0,sqrt0,sqrt0/2,sqrt0/50]
titles=[rf'$2k{sign}\Delta\Omega$' for sign in ['\gg','>','=','<','\ll']]

for a,k,title in zip(ax,k0,titles):
    L.kex=sl.Tools.twoSite_kex(tc=1/(2*k),p1=.3)
    OneD(L,ax=a)
    a.set_title(title)
_=fig.tight_layout()


# We repeat the eigenvalue calculation.

# In[18]:


rho=sl.Rho('S0x','S0p')
rho_r,_=rho.ReducedSetup(L.U(.1))
L_r=rho_r.L   # The reduced Liouvillian is not returned directly
for k in k0:
    L.kex=sl.Tools.twoSite_kex(tc=1/(2*k),p1=.3)
    Lambda=np.linalg.eig(L_r[0].L(0))[0]
    print(f'k={k:.1} s^-1,'+\
          f'T2={-1e3/Lambda[0].real:5.2f} ms, T2={-1e3/Lambda[1].real:5.2f} ms,'+\
          f'f={-Lambda[0].imag/2/np.pi:5.1f} Hz, f={-Lambda[1].imag/2/np.pi:5.1f} Hz')


# ### Exercise 1.4
# What do you see that's different about the results of the eigenvalue calculation above vs. symmetric exchange.
# ```{toggle}
# The critical point is that there are always two unique frequencies for asymmetric exchange. In symmetric exchange, when $2k>\Delta\Omega$, a coalescence condition occured and the two peaks both had the same frequency. In asymmetric exchange, no coalescence exists. One peak just becomes increasingly weak until it is effectively invisible.
# ```

# In chemical exchange in 1D, in principle, we can extract populations of the two states and rates of exchange between them, simply by fitting the spectrum. However, other sources of relaxation (e.g. tumbling in solution) will contribute to the broadening, complicating the extraction of exchange rates and populations from the spectra alone. In this case, we usually need a reference sample which should have the same $T_2$ without exchange present. This is not always available, so one may need to consider other approaches.
# 
# A second challenge is what to do when there is not significant broadening in the 1D experiment; indeed, unless $2k\approx\Delta\Omega$, this is the case.

# ## Chemical Exchange in 2D

# In the oneD spectra above, we get significant broadening for the case that rate of exchange ($2k$) is on the order of the size in chemical shift. However, what should we do if this is not the case? The first situation is if motion is much slower. Then, we obtain narrow lines, where broadening is likely dominated by other sources of relaxation so that we cannot easily extract exchange rates.
# 
# In this case, we may use Exchange Spectroscopy (EXSY). This is a 2D experiment, where we evolve in the indirect dimension, flip the magnetization to the *z*-axis, wait to allow exchange to occur, and then flip magnetization back to the *x*-axis to observe the emergence of new frequencies.
# 
# We include some $T_2$ relaxation to broaden the spectra slightly.

# In[19]:


L.kex=sl.Tools.twoSite_kex(tc=.2,p1=0.6)
L.add_relax('T2',i=0,T2=.05)

_=EXSY(L,dly=1,buildup=False) #dly is the length of time between flip up and flip down


# ### Exercise 2.1
# Write formulas for the peak integrals for a delay much longer than the length of the correlation time of exchange. That is, what is $I_{1,1}$, $I_{1,2}$, $I_{2,1}$, $I_{2,2}$ as a function of $p_1$ and $p_2$, where $I_{1,2}$ is the peak corresponding to magnetization starting in state one and ending in state 2, for example, and $p_1$ and $p_2$
# 
# 
# We can take the populations as probabilities, and so the probability of starting in state 1 and ending in state one is simply $p_1*p_1$. By analogy:
# 
# ```{toggle}
# $$
# \begin{eqnarray}
# I_{1,1}&=&p_1^2 \\ 
# I_{1,2}&=&p_1p_2 \\
# I_{2,1}&=&p_1p_2 \\
# I_{2,2}&=&p_2^2
# \end{eqnarray}
# $$
# ```

# ### Exercise 2.2
# 
# For a given delay, $\tau$, and exchange rate constants $k_{12}$ and $k_{21}$, give the peak intensities.
# 
# Hint: Give your answer in terms of $p_1$, $p_2$, and $k$ ($k=k_1+k_2$).
# 
# 
# $$
# \begin{eqnarray}
# I_{1,1}(\tau)&=&p_1^2+(p_1-p_1^2)\exp(-k\tau) \\
# I_{1,2}(\tau)&=&p_1p_2(1-\exp(-k\tau)) \\
# I_{2,1}(\tau)&=&p_1p_2(1-\exp(-k\tau)) \\
# I_{2,2}(\tau)&=&p_2^2+(p_2-p_2^2)\exp(-k\tau)
# \end{eqnarray}
# $$

# Below, we plot the buildup of peak intensities for our system, where we observe the expected behavior. Dotted lines show the result from the formulas above.

# In[20]:


dly=1
ax=EXSY(L,dly=dly)
t=np.linspace(0,dly,20)
k=L.kex[0,1]+L.kex[1,0]
p1=1/(1+L.kex[1,0]/L.kex[0,1])
p2=1-p1
ax.plot(t*1e3,p1**2+(p1-p1**2)*np.exp(-t*k),linestyle=':',color='grey')
ax.plot(t*1e3,p2**2+(p2-p2**2)*np.exp(-t*k),linestyle=':',color='grey')
ax.plot(t*1e3,p1*p2*(1-np.exp(-t*k)),linestyle=':',color='grey')


# The agreement is not quite perfect. Indeed, some exchange may occur during the indirect and direct dimensions, so that there is some cross-peak amplitude even when the delays are zero. Nonetheless, this is a good model of the buildup.

# ### Exercise 2.3
# When, roughly, is exchange too fast (or two slow) for the EXSY experiment? What are the problems that arise? Edit `k` in the cell below. When `k` is large, you may want to shorten the total delay in the second call to EXSY. Note that we have added $T_1$ relaxation to the Liouvillian (why is this relevant? Can you edit the buildup formulas accordingly?)

# In[21]:


k=.15
L.clear_relax()
L.add_relax('T2',i=0,T2=.05)
L.add_relax('T1',i=0,T1=2)
L.kex=sl.Tools.twoSite_kex(tc=1/(2*k),p1=0.6)

fig=plt.figure(figsize=[10,4])
ax=[fig.add_subplot(1,2,1,projection='3d'),fig.add_subplot(1,2,2)]
EXSY(L,dly=1,buildup=False,ax=ax[0])
_=EXSY(L,dly=5,buildup=True,ax=ax[1])


# ```{toggle}
# Around $k>1500$, the 2D spectrum becomes too broad to separate the individual peaks and so we can no longer easily extract exchange rates from the EXSY buildup curves. On the other hand, if $k<0.15$, the peaks decay faster than any significant transfer occurs. Then, the formulas for the peak intensities are
# 
# $$
# \begin{eqnarray}
# I_{1,1}(\tau)&=&(p_1^2+(p_1-p_1^2)\exp(-k\tau))\cdot\exp(-\tau/T_1) \\
# I_{1,2}(\tau)&=&p_1p_2(1-\exp(-k\tau))\cdot\exp(-\tau/T_1) \\
# I_{2,1}(\tau)&=&p_1p_2(1-\exp(-k\tau))\cdot\exp(-\tau/T_1) \\
# I_{2,2}(\tau)&=(&p_2^2+(p_2-p_2^2)\exp(-k\tau))\cdot\exp(-\tau/T_1)
# \end{eqnarray}
# $$
# 
# Then, the relevant parameters are $\Delta\Omega$, where if $k>>\DeltaOmega$, peaks become too broad, and if $k<<1/T_1$, then signal intensity is lost too quickly to use EXSY
# 
# ```

# In[22]:


fig=plt.figure(figsize=[12,8])
ax=[fig.add_subplot(2,2,1,projection='3d'),fig.add_subplot(2,2,2),
   fig.add_subplot(2,2,3,projection='3d'),fig.add_subplot(2,2,4)]

k=1500
L.kex=sl.Tools.twoSite_kex(tc=1/(2*k),p1=0.6)
EXSY(L,dly=1,buildup=False,ax=ax[0])
_=EXSY(L,dly=.005,buildup=True,ax=ax[1])

k=0.15
L.kex=sl.Tools.twoSite_kex(tc=1/(2*k),p1=0.6)
EXSY(L,dly=1,buildup=False,ax=ax[2])
_=EXSY(L,dly=5,buildup=True,ax=ax[3])


# ## CEST
# Our experiments up until now assume a relatively large population on both states in exchange. However, what can we do if one population is too small to observe as a resonance in the spectrum, we may use the Chemical Exchange Saturation Transfer experiment, or CEST. This experiment applies a saturating field onto longitudinal magnetization followed by detection. If the saturating field is on-resonance with the main peak, then as expected, the signal is reduced. However, if a peak corresponding to a low population is on-resonance, then the main peak also saturates. This is demonstrated below, where we first plot the 1D experiment, followed by the CEST saturation profile.

# In[23]:


L.clear_relax()
L.add_relax('T2',i=0,T2=.05)
L.add_relax('T1',i=0,T1=2)
L.add_relax('recovery') #Causes magnetization to recover to thermal equilibrium
L.kex=sl.Tools.twoSite_kex(tc=1e-2,p1=0.98)

fig,ax=plt.subplots(1,2,figsize=[9,4])
OneD(L,ax=ax[0])
_=CEST(L,dly=1,v1=25,ax=ax[1])


# The above simulation has set the second peak to be just 2% of the total signal. However, the peak is easily visible in the saturation profile. 

# ### Exercise 4.1
# Below, `dly` gives the length of the CEST saturation period in seconds, and `v1` gives the strength of the saturation field in Hz. Fix the experiment! What should one consider when setting these parameters?

# In[24]:


_=CEST(L,dly=.01,v1=500)


# ```{toggle}
# Set `dly`=1 s, and `v1`=25 Hz. 
# 
# Using too much power saturates the whole spectrum, and so we lose detailed features. However, we need to apply the field long enough that saturation actually occurs. 
# 
# ```

# ### Exercise 4.2
# How does the CEST experiment depend on the rate of exchange? Adjust `k` below to determine where CEST works or fails. Don't forget to fix `dly` and `v1` based on your investigation above.

# In[25]:


k=1
L.kex=sl.Tools.twoSite_kex(tc=1/(2*k),p1=0.98)

_=CEST(L,dly=5,v1=25)


# ```{toggle}
# When `k` approches ~2000 s$^{-1}$, peaks begin to coalesce, making it more difficult to identify location of the second peak. On the other hand, if `k` becomes too small, then saturation on the weak resonance does not get transferred to the main resonance quickly enough, occuring when `k` approaches 0.5 s$^{-1}$.
# ```

# ## Matching Game
# Below, we have prepared four Liouvillian matrices (objects), corresponding to different populations and correlation times. Which experiment would you apply for each Liouvillian? Try out each experiment (use CPMG or $R_{1p}$, but not both, since their performance is very similar)
# 
# Hint: There is some overlap for which experiment could work with each Liouvillian, but there is a "best" one-to-one assigment. You can adjust the input parameters for the experiments, but this should not be required.

# In[102]:


L1,L2,L3,L4=[sl.Liouvillian(ex0,ex1) for _ in range(4)]
for L,(k,p1) in zip([L1,L2,L3,L4],pars):
    L.add_relax('T2',i=0,T2=0.05)
    L.add_relax('T1',i=0,T1=2)
    L.kex=sl.Tools.twoSite_kex(k,p1)


# Reveal the cell below for our suggested best answer

# In[133]:


fig,ax=plt.subplots(2,2,figsize=[9,9])
ax=ax.flatten()
for k,(L,exp,a) in enumerate(zip([L1,L2,L3,L4],[CEST,R1p,EXSY,OneD],ax)):
    exp(L,ax=a)
    a.set_title(fr'L{k+1}$\rightarrow${exp.__name__}')
fig.tight_layout()

for k,(L,(tc,p1)) in enumerate(zip([L1,L2,L3,L4],pars)):
    print(f'L{k+1}: k={1/(2*tc):.2}, p1={p1:4.2f}')

