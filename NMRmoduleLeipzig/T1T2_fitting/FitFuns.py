import numpy as np                #Math tools
import matplotlib.pyplot as plt   #Plotting tools
from scipy.optimize import least_squares as lsq  #Non-linear least squares

def ExpFit(t,I,mode='IR'):
    """
    Fits data, specified by a time axis and an intensity, to an exponential decay.

    Modes are:
    'IR' (inversion recovery)
        I=A(1-2*exp(-t/tau)
    'SR' (saturation recovery)
        I=A(1-exp(-t/tau)
    'T2' (T2 decay)
        I=A*exp(-t/tau)

    returns
    A : Amplitude of signal
    tau : Time constant of decay
    model : Functional form of model
    """
    t=np.array(t)
    I=np.array(I)
    if mode=='IR':
        def model(X,t):
            return X[0]-2*X[0]*np.exp(-t/X[1])
    elif mode=='SR':
        def model(X,t):
            return X[0]*(1-np.exp(-t/X[1]))
    else:
        def model(X,t):
            return X[0]*np.exp(-t/X[1])

    def fun(X,t,I):
        return model(X,t)-I
        
    
    if 'IR':
        b=np.argmin(np.abs(I))
        X0=(-I[0],t[b])
    elif 'SR':
        b=np.argmin(np.abs(I/I[-1]-.63))
        X0=(I[-1],t[b])
    else:
        b=np.argmin(np.abs(I/I[0]-.37))
        X0=(I[0],t[b])
    out=lsq(fun,X0,args=(t,I))
    return *out['x'],model

def PlotFit(t,I,mode='IR'):

    assert len(I)==len(t),f"The number of elements in t ({len(t)}) must be the same as the number of elements in I ({len(I)})"
    
    A,tau,model=ExpFit(t,I,mode)
    ax=plt.subplots()[1]
    t0=np.linspace(t[0],t[-1],100)
    Ifit=model([A,tau],t0)
    ax.plot(t0,Ifit,label='fit',color='black')
    ax.scatter(t,I,label='input',color='red')
    ax.legend()
    ax.set_xlabel('t / s')
    ax.set_ylabel('Intensity / a.u.')
    return ax,tau


def T1fit(t,I):
    """
    Extracts the T1 from a saturation recovery experiment, assuming the following function form
    
    I=A(1-2*exp(-t/T1)

    Arguments:
    t : Time points acquired (in seconds)
    I : Amplitude of signal at each time point

    returns T1 (in seconds)
    
    """
    ax,T1=PlotFit(t,I,mode='IR')
    ax.text(np.mean(ax.get_xlim()),np.mean(ax.get_ylim()),fr'$T_1$ = {T1:.2f} s')
    return T1

def T2fit(t,I):
    """
    Extracts the T2 from a CPMG experiment, assuming the following function form
    
    I=A*exp(-t/T2)

    Arguments:
    t : Time points acquired (in seconds)
    I : Amplitude of signal at each time point

    returns T2 (in seconds)
    
    """
    
    ax,T2=PlotFit(t,I,mode='T2')
    ax.text(np.mean(ax.get_xlim()),np.mean(ax.get_ylim()),fr'$T_2$ = {T2:.2f} s')
    return T2