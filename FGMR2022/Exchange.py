#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 12:47:19 2022

@author: albertsmith
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 11:00:12 2022

@author: albertsmith
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cv2



class Exchange():
    def __init__(self,v1:float=1,v2:float=10,k12:float=5,k21:float=None):
        """
        Initializes the Exchange class, with frequencies of the two states and
        corresponding exchange rates

        Parameters
        ----------
        v1 : float, optional
            Chemical shift (Hz) of state 1. The default is 2.
        v2 : float, optional
            Chemical shift (Hz) of state 2. The default is 15.
        k12 : float, optional
            Rate constant for exchange from state 1 to 2. The default is 2.
        k21 : float, optional
            Rate constant for exchange from state 2 to 1. The default is None,
            in which case the value will be set to k12.
            
        Returns
        -------
        None.

        """
        self.v1=v1
        self.v2=v2
        self.k12=k12
        self._k21=k21
        self._DV=None
        self._frind=0
        self._moviefilename=None
        self.fig=None
        self.rate_scaling=10
        self.fr=15
    
    @property
    def k21(self):
        return self.k12 if self._k21 is None else self._k21
    
    def __setattr__(self,name,value):
        if name in ['v1','v2','k12','k21']:
            self._DV=None
        if name=='k21':
            self._k21=value
            return
        super().__setattr__(name,value)
    
    @property
    def exmat(self) -> np.array:
        """
        Returns the exchange matrix for exchange of the two states.

        Returns
        -------
        np.array
            Exchange matrix.

        """
        return np.array([[-self.k12,self.k21],[self.k12,-self.k21]])
    
    @property
    def LiouvilleMat(self) -> np.array:
        """
        Returns exchange matrix for one coherence of the full Liouville matrix,
        i.e. just the exchange matrix plus the frequencies along the diagonal.

        Returns
        -------
        np.array
            Exchange matrix + coherence frequencies

        """
        return self.exmat+np.diag([1j*self.v1*2*np.pi,1j*self.v2*2*np.pi])
    
    @property
    def pop(self) -> np.array:
        """
        Populations of state 1 and state 2

        Returns
        -------
        pop : np.array
            Populations of state 1 and state 2.

        """
        D,V=np.linalg.eig(self.exmat)
        i=np.argsort(np.abs(D))
        pop=V[:,i[0]].real
        pop/=pop.sum()
        return pop
    
    @property
    def k(self) -> float:
        """
        Exchange rate constant, i.e. the non-zero eigenvalue of the exchange
        matrix

        Returns
        -------
        float
            Exchange rate constant.

        """
        D,_=np.linalg.eig(self.exmat)
        i=np.argsort(np.abs(D))
        return -D[i[1]].real
    
    @property
    def kcoal(self) -> float:
        """
        Returns the rate constant of exchange (symmetric) that results in
        coalescence of the signal
        
        k=pi*|v1-v2|

        Returns
        -------
        float
            Rate constant yielding coalescence

        """
        return np.pi*np.abs(self.v1-self.v2)
    
    @property
    def tc(self) -> float:
        """
        Correlation time, i.e. the inverse of the non-zero eigenvalue of the 
        exchange matrix.

        Returns
        -------
        float
            Correlation time of exchange matrix
        """
        return 1/self.k
    
    @property
    def LMeig(self) -> tuple:
        """
        Yields the eigenvectors and eigenvalues of the Liouville matrix

        Returns
        -------
        tuple
            DESCRIPTION.

        """
        if self._DV is None:
            self._DV=np.linalg.eig(self.LiouvilleMat)
        return self._DV
    
    @property
    def v(self) -> np.array:
        """
        Returns the observed frequencies for the two signals (may be identical)

        Returns
        -------
        np.array
            Frequencies of the two signals

        """
        return self.LMeig[0].imag/2/np.pi
    
    @property
    def T2(self) -> np.array:
        """
        Returns the relaxation rate constants for the two signals

        Returns
        -------
        np.array
            Relaxation rate constants

        """
        
        return 1/self.FWHM/np.pi
    
    @property
    def FWHM(self) -> np.array:
        """
        Returns the linewidths of the two signals (Full width at half maximum)

        Returns
        -------
        np.array
            Full width at half maximum

        """

        return -self.LMeig[0].real/np.pi
    
    
    @property
    def A(self) -> np.array:
        """
        Returns the amplitudes corresponding to the two signals. Note this is
        the time-domain amplitude, corresponding to the integral in frequency
        domain

        Returns
        -------
        np.array
            Amplitude

        """
        # A=self.LMeig[1].T@self.pop
        # A/=A.sum().real #Is this correct?
        # return A
        
        
        
        V=self.LMeig[1]
        
        return V.sum(0)*np.linalg.solve(V,self.pop)
        
    
    def t(self,tf:float=None,nt:int=None):
        """
        Generate a time axis, with optional default values

        Parameters
        ----------
        tf : float, optional
            Total length of the trajectory. The default is None, which will
            complete 10 cycles of the slower frequency.
        nt : int, optional
            Total number of time points. The default is None, which will place
            50 time points within one cycle of the faster frequency.

        Returns
        -------
        None.

        """
        if tf is None:
            tf=1/np.max(self.v)*8
        
        if nt is None:
            nt=int(tf*np.max(self.v)*60)+1

        return np.linspace(0,tf,nt)
    
    def f(self,tf:float=None,nt:int=None):
        """
        Returns the frequency axis for a given timescale axis, assuming zero-filling
        to 2x the length of the trajectory

        Parameters
        ----------
        tf : float, optional
            Total length of the trajectory. The default is None, which will
            complete 10 cycles of the slower frequency.
        nt : int, optional
            Total number of time points. The default is None, which will place
            50 time points within one cycle of the faster frequency.

        Returns
        -------
        np.array
            frequency axis

        """
        t=self.t(tf=tf,nt=nt)
        dt=t[1]-t[0]
        f=np.linspace(-1,1,len(t)*2)/dt/2
        f+=(f[0]-f[1])/2
        return f
        
    
    def calc_traj(self,state0:int=None,tf:float=None,nt:int=None):
        """
        Calculates the trajectory and spectrum for an ensemble of spins
        

        Parameters
        ----------
        tf : float, optional
            Total length of the trajectory. The default is None, which will
            complete 10 cycles of the slower frequency.
        nt : int, optional
            Total number of time points. The default is None, which will place
            50 time points within one cycle of the faster frequency.

        Returns
        -------
        None.

        """
        t=self.t(tf=tf,nt=nt)
        
        if state0 is None:
            I=(self.A*np.exp(np.atleast_2d(t).T@np.atleast_2d(2*np.pi*1j*self.v-1/self.T2))).sum(1)
        else:
            state0-=1
            I=(self.A[state0]*np.exp(np.atleast_2d(t).T@np.atleast_2d(2*np.pi*1j*self.v[state0]-1/self.T2[state0]))).sum(1)
            
        S=np.fft.fftshift(np.fft.fft(np.concatenate(([I[0]/2],I[1:])),n=len(I)*2))
        return {'t':t,'I':I,'f':self.f(tf=tf,nt=nt),'S':S}
    
    def single_spin_traj(self,state0:int=1,tf:float=None,nt:int=None):
        """
        Calculates the trajectory of a single spin hopping between the two
        states, with randomly calculated trajectory. Starts in state 1 or 2,
        depending on state0 setting.

        Parameters
        ----------
        state0 : int, optional
            Initial state of the spin (1 or 2). The default is 1 
        tf : float, optional
            Total length of the trajectory. The default is None, which will
            complete 10 cycles of the slower frequency.
        nt : int, optional
            Total number of time points. The default is None, which will place
            50 time points within one cycle of the faster frequency.

        Returns
        -------
        dict
            t : time axis
            state : state vector (1 or 2) for each time point
            I : Complex  signal (time domain)
            f : frequency axis
            S : Frequency domain signal

        """
        
        
        t=self.t(tf=tf,nt=nt)
        nt=len(t)

        dt=t[1]
        p1,p2=1-np.exp(-dt*np.array([self.k12,self.k21]))
        hop1=np.random.rand(nt)<p1
        hop2=np.random.rand(nt)<p2
        
        pos=0
        state=np.zeros(nt,dtype='int8')
        state[0]=state0
        while True:
            if state[pos]==1:
                if not(any(hop1[pos+1:])):
                    state[pos:]=1
                    break
                i=np.argwhere(hop1[pos+1:])[0,0]+1
                state[pos:pos+i]=1
                state[pos+i]=2
                pos+=i
            else:
                if not(any(hop2[pos+1:])):
                    state[pos:]=2
                    break
                i=np.argwhere(hop2[pos+1:])[0,0]+1
                state[pos:pos+i]=2
                state[pos+i]=1
                pos+=i
        
        v=np.zeros(nt)
        v[state==1]=self.v1
        v[state==2]=self.v2
        v=np.concatenate(([0],v[:-1]))
        I=np.exp(1j*2*np.pi*np.cumsum(v*dt))
                
        S=np.fft.fftshift(np.fft.fft(np.concatenate(([I[0]/2],I[1:])),n=len(I)*2))
        return {'t':t,'state':state,'I':I,'f':self.f(tf=tf,nt=nt),'S':S}
    
    def accumulate_1spin(self,state0:int=None,n:int=100,tf:float=None,nt:int=None):
        """
        Accumulates 1 spin spectra with random trajectories. One may specify
        the starting state (state0), or leave to None, which will then accumulate
        both starting states, proportionally to their population

        Parameters
        ----------
        state0 : int, optional
            Initial state of the spin (1 or 2). The default is 1 
        n : int, optional
            Number of 1 spin spectra to accumulate. The default is 100.
        tf : float, optional
            Total length of the trajectory. The default is None, which will
            complete 10 cycles of the slower frequency.
        nt : int, optional
            Total number of time points. The default is None, which will place
            50 time points within one cycle of the faster frequency.

        Returns
        -------
        None.

        """
        t=self.t(tf=tf,nt=nt)
        out={'t':t,'I':np.zeros([n,len(t)]),'Iacc':None,
             'f':self.f(tf=tf,nt=nt),'S':np.zeros([n,2*len(t)]),'Sacc':None}
        
        for k in range(n):
            s0=(1 if np.random.rand(1)<self.pop[0] else 2) if state0 is None else state0
            temp=self.single_spin_traj(s0,tf=tf,nt=nt)
            out['I'][k]=temp['I']
            out['S'][k]=temp['S']
        out['Iacc']=np.cumsum(out['I'],axis=0)
        out['Sacc']=np.cumsum(out['S'],axis=0)
        return out
    

            
#%% Movies  
    def single_spin_movie(self,state0:int=1,tf:float=None,nt:int=None,ax=None):
        """
        Creates of movie of the trajectory of a single spin

        Parameters
        ----------
        state0 : int, optional
            Initial state of the spin (1 or 2). The default is 1 
        tf : float, optional
            Total length of the trajectory. The default is None, which will
            complete 10 cycles of the slower frequency.
        nt : int, optional
            Total number of time points. The default is None, which will place
            50 time points within one cycle of the faster frequency.
        ax : plt.axis, optional
            axis on which we plot the trajectory. Generated if not provided.
            The default is None

        Returns
        -------
        None.

        """
        if ax is None:
            self.fig=plt.figure()
            ax=self.fig.add_subplot(111)
        
        ax.set_xlim([-1.1,1.1]) 
        ax.set_ylim([-1.1,1.1])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal','box')
        hdl=ax.plot([0,1],[0,0],color='black')[0]
        
        traj=self.single_spin_traj(state0=state0,tf=tf,nt=nt)
        
        hdl0=list()
        hdl1=list()
        for I0,I1,state in zip(traj['I'][:-1],traj['I'][1:],traj['state']):
            self.single_frame()
            hdl.set_xdata([0,I1.real])
            hdl.set_ydata([0,I1.imag])
            hdl0.append(ax.plot([I0.real,I1.real],[I0.imag,I1.imag],color='blue' if state==1 else 'red')[0])
            if I0.imag<0 and I1 .imag>=0:
                for h in hdl1:h.remove()
                for h in hdl0:h.set_alpha(.25)
                hdl1=hdl0
                hdl0=list()
                
        return traj
    
    def single_spin_v_rate(self,state0:int=1,tf:float=None,nt:int=None):
        """
        Creates a movie with 3 single-spin trajectories, where the rate of motion
        is slow (k=1/5*k_coalescence), intermediate (k=k_coalescence), and fast
        (k=5*k_coalescence).

        Parameters
        ----------
        state0 : int, optional
            Initial state of the spin (1 or 2). The default is 1 
        tf : float, optional
            Total length of the trajectory. The default is None, which will
            complete 10 cycles of the slower frequency.
        nt : int, optional
            Total number of time points. The default is None, which will place
            50 time points within one cycle of the faster frequency.

        Returns
        -------
        None.

        """
        
        self.fig,ax=plt.subplots(3,2)
        rsc=self.rate_scaling
        titles=[r'Slow motion ($k=$'+'{:.1f}'.format(1/rsc)+r'$*k_{coalescence})$',r'Coalescence ($k=k_{coalescence})$',
                r'Fast motion ($k=$'+'{:.0f}'.format(rsc)+r'$*k_{coalescence})$']
        # hdl=list()
        for a,title in zip(ax,titles):
            a[0].set_xlim([-1.1,1.1]) 
            a[0].set_ylim([-1.1,1.1])
            a[0].set_xlabel('x')
            a[0].set_ylabel('y')
            a[0].set_aspect('equal','box')
            a[1].set_title(title)
            width=6*np.abs(self.v2-self.v1)
            a[1].set_xlim(self.v.mean()+np.array([-width/2,width/2]))
            a[1].set_ylim([-1.1,1.1])
            a[1].set_xlabel(r'$\delta$ / Hz')
            # hdl.append(a[0].plot([0,1],[0,0],color='black')[0])
        self.fig.set_size_inches([7.62, 7.01])
        self.fig.tight_layout()
        
        self.k21=None
        for k,a in zip(np.array([1/rsc,1,rsc])*self.kcoal,ax):
            self.k12=k
            # traj=self.single_spin_movie(state0=state0,tf=tf,nt=nt,ax=a[0])
            
            #Insert new code here
            a[0].set_xlim([-1.1,1.1]) 
            a[0].set_ylim([-1.1,1.1])
            a[0].set_xlabel('x')
            a[0].set_ylabel('y')
            a[0].set_aspect('equal','box')
            hdl=a[0].plot([0,1],[0,0],color='black')[0]

            
            traj=self.single_spin_traj(state0=state0,tf=tf,nt=nt)
            S=np.fft.fftshift(np.fft.fft(np.concatenate(([traj['I'][0]/2],traj['I'][1:2])),n=len(traj['f'])))
            hdls=a[1].plot(traj['f'],S.real/S.real.max())[0]
            
            hdl0=list()
            hdl1=list()
            for m,(I0,I1,state) in enumerate(zip(traj['I'][:-1],traj['I'][1:],traj['state'])):
                self.single_frame()
                hdl.set_xdata([0,I1.real])
                hdl.set_ydata([0,I1.imag])
                hdl0.append(a[0].plot([I0.real,I1.real],[I0.imag,I1.imag],color='blue' if state==1 else 'red')[0])
                if I0.imag<0 and I1 .imag>=0:
                    for h in hdl1:h.remove()
                    for h in hdl0:h.set_alpha(.25)
                    hdl1=hdl0
                    hdl0=list()
                   
                S=np.fft.fftshift(np.fft.fft(np.concatenate(([traj['I'][0]/2],traj['I'][1:m+1])),n=len(traj['f'])))
                hdls.set_ydata(S.real/S.real.max())
            #To here
            
            # a[1].plot(traj['f'],traj['S'].real/traj['S'].real.max())  #Comment here
            for _ in range(15):self.single_frame()
            
        
            # if len(hdl0)>30:hdl0[len(hdl0)-30].set_alpha(.5)
            # if len(hdl0)>60:hdl0[len(hdl0)-60].remove()
            
    def movie_accum(self,n:int=300,step:int=2,state0:int=None,tf:float=None,nt:int=None,ax=None):
        """
        

        Parameters
        ----------
        n : int
        state0 : int, optional
            Initial state of the spin (1 or 2). The default is 1 
        tf : float, optional
            Total length of the trajectory. The default is None, which will
            complete 10 cycles of the slower frequency.
        nt : int, optional
            Total number of time points. The default is None, which will place
            50 time points within one cycle of the faster frequency.

        Returns
        -------
        None.

        """
        
        traj=self.accumulate_1spin(state0=state0,n=n,tf=tf,nt=nt)
        
        if ax is None:
            self.fig=plt.figure()
            ax=self.fig.add_subplot(111)
        
        out=self.calc_traj(state0=state0,tf=tf,nt=nt)
        S=out['S']/out['I'].shape[0]
        sc=1/S.max()
        sc1=sc if self.k12<=self.kcoal else sc*out['S'].max()/traj['Sacc'][-1].max()*n
        
        hdl=ax.plot(traj['f'],traj['Sacc'][0]/traj['I'].shape[1]*sc1,color='black')[0]
        width=6*np.abs(self.v2-self.v1)
        ax.set_xlim(self.v.mean()+np.array([-width/2,width/2]))
        ax.set_ylim([-.5,2])
        ax.set_xlabel(r'$\delta$ / Hz')
        
        for k in range(0,n,step):
            self.single_frame()
            hdl.set_ydata(traj['Sacc'][k]/(k+1)/traj['I'].shape[1]*sc1)
        self.single_frame()
        # out=self.calc_traj(state0=state0,tf=tf,nt=nt)
        # S=out['S']/out['I'].shape[0]
        ax.plot(traj['f'],S*sc,color='red')
        self.single_frame()
        
            

    def movie_accum_v_rate(self,n:int=300,step:int=2,state0:int=None,tf:float=None,nt:int=None):
        """
        

        Parameters
        ----------
        n : int
        state0 : int, optional
            Initial state of the spin (1 or 2). The default is 1 
        tf : float, optional
            Total length of the trajectory. The default is None, which will
            complete 10 cycles of the slower frequency.
        nt : int, optional
            Total number of time points. The default is None, which will place
            50 time points within one cycle of the faster frequency.

        Returns
        -------
        None.

        """
        
        self.fig,ax=plt.subplots(1,3)
        rsc=self.rate_scaling
        titles=[r'Slow motion ($k=$'+'{:.1f}'.format(1/rsc)+r'$*k_{coalescence})$',r'Coalescence ($k=k_{coalescence})$',
                r'Fast motion ($k=$'+'{:.0f}'.format(rsc)+r'$*k_{coalescence})$']
        
        for a,title in zip(ax,titles):
            a.set_title(title)
            width=6*np.abs(self.v2-self.v1)
            a.set_xlim(self.v.mean()+np.array([-width/2,width/2]))
            a.set_ylim([-1.1,1.1])
            a.set_xlabel(r'$\delta$ / Hz')
            
        self.fig.set_size_inches([10.5,3.9])
        self.fig.tight_layout()
        
        self.k21=None
        for k,a in zip(np.array([1/rsc,1,rsc])*self.kcoal,ax):
            self.k12=k
            self.movie_accum(n=n,step=step,state0=state0,tf=tf,nt=nt,ax=a)
            
    def movie_ksweep(self,ratio:list=np.logspace(-1.5,1.5,75),tf=None,nt=None,ax=None):
        """
        Creates a movies sweeping through the coalescence condition. Input is a
        list of ratios to use relative to the coalescence frequency

        Parameters
        ----------
        ratio : list, optional
            DESCRIPTION. The default is np.logspace(-1.5,1.5,75).

        Returns
        -------
        None.

        """
        
        if ax is None:
            self.fig=plt.figure()
            ax=self.fig.add_subplot(111)
        ratio=np.array(ratio)
        
        self.k12=self.kcoal*ratio[0]
        if nt is None:nt=self.t().__len__()
        if tf is None:tf=self.t()[-1]
        out=self.calc_traj(tf=tf,nt=nt)
        hdl=ax.plot(out['f'],out['S'],color='red',linewidth=2)[0]
        ax.set_xlim(4*np.abs(self.v1-self.v2)*np.array([-.5,.5])+(self.v1+self.v2)/2)
        ax.set_xlabel(r'$\delta$ / Hz')
        ax.set_title(r'$k$ = '+'{0:.1f}'.format(self.k12)+r' s$^{-1}$, $\Delta\Omega/2\pi$'+' = {0:.1f} Hz'.format(np.abs(self.v1-self.v2)))
        self.k12=self.kcoal*ratio[-1]
        ymax=self.calc_traj(tf=tf,nt=nt)['S'].max()
        ax.set_ylim([-0.1*ymax,1.05*ymax])
        
        self.single_frame()
        for rat in ratio[1:]:
            self.k12=self.kcoal*rat
            out=self.calc_traj(tf=tf,nt=nt)
            hdl.set_ydata(out['S'])
            ax.set_title(r'$k$ = '+'{0:.1f}'.format(self.k12)+r' s$^{-1}$, $\Delta\Omega/2\pi$'+' = {0:.1f} Hz'.format(np.abs(self.v1-self.v2)))
            self.single_frame()
        
        
    #%% Multi-spin visualizaton
    def make3Daxis(self,ax=None):
        """
        Makes a 3D axis for plotting

        Parameters
        ----------
        ax : TYPE, optional
            Axis to edit for the 3D axis. The default is None.

        Returns
        -------
        None.

        """
        if ax is None:ax = plt.figure().add_subplot(111, projection='3d')
        
        self.fig=ax.figure
        t=np.linspace(0,2*np.pi,150)
        
        # ax.plot3D(np.zeros(t.shape),np.sin(t),np.cos(t),color='grey')
        
        
        # ax.plot3D(np.sin(t),np.zeros(t.shape),np.cos(t),color='grey')
        
        rat=ax.get_box_aspect()[1]/ax.get_box_aspect()[2]
        ax.set_xlim([-rat,rat])
        ax.set_ylim([-rat,rat])
        ax.set_zlim([-1,1])
        ax.plot3D([0,1],[0,0],[0,0],color='red',linewidth=2)
        ax.plot3D([0,-1],[0,0],[0,0],color='red',linewidth=2,alpha=.25)
        ax.plot3D([0,0],[0,1],[0,0],color='green',linewidth=2)
        ax.plot3D([0,0],[0,-1],[0,0],color='green',linewidth=2,alpha=.25)
        ax.plot3D([0,0],[0,0],[0,1],color='blue',linewidth=2)
        ax.plot3D([0,0],[0,0],[0,-1],color='blue',linewidth=2,alpha=.25)
        ax.view_init(elev=31.5,azim=46.7)

        
        for k in np.linspace(-1,1,5):
            ax.plot3D(np.cos(t)*np.sqrt(1-k**2),np.sin(t)*np.sqrt(1-k**2),np.ones(t.shape)*k,color='grey')
        for k in np.linspace(0,np.pi,7):
            ax.plot3D(np.sin(t)*np.cos(k),np.sin(t)*np.sin(k),np.cos(t),color='grey')
        ax.plot3D(np.cos(t),np.sin(t),np.zeros(t.shape),color='black')
        return ax    
        
    def movie_hahn_static(self,v:list=None,ax=None):
        """
        Generates trajectories with different frequencies and plots these, 
        including a pi pulse to refocus them

        Parameters
        ----------
        v : TYPE, optional
            DESCRIPTION. The default is np.linspace(-25,25,51).

        Returns
        -------
        None.

        """
        
        if ax is None:
            self.fig=plt.figure()
            ax=self.fig.add_subplot(1,1,1,projection='3d')
        
        self.make3Daxis(ax)
        # v=np.linspace(self.v1,self.v2,11) if v is None else np.array(v)
        v=np.random.rand(11)*(self.v2-self.v1)+self.v1

        
        t=np.atleast_2d(np.linspace(0,1/v.mean(),225)).T
        v=np.atleast_2d(v)
        
        M=np.zeros([3,len(v),len(t)*2+10])
        
        M0=np.zeros([3,len(v)])
        M0[0]=1
        c,s=np.cos(2*np.pi*t@v),np.sin(2*np.pi*t@v)
        M=np.array([c*M0[0]-s*M0[1],s*M0[0]+c*M0[1],M0[2]*np.ones(c.shape)])
        M0=[np.atleast_2d(x) for x in M[:,-1]]
        c1,s1=[np.atleast_2d(getattr(np,x)(np.linspace(0,np.pi,10))).T for x in ['cos','sin']]
        M1=np.array([np.ones(c1.shape)@M0[0],c1@M0[1],-s1@M0[1]])       
        M0=M1[:,-1]
        M2=np.array([c*M0[0]-s*M0[1],s*M0[0]+c*M0[1],M0[2]*np.ones(c.shape)])
        M=np.moveaxis(np.concatenate((M,M1,M2),axis=1),0,2)
        
        cmap=plt.get_cmap('tab10')
        hdl=[ax.plot3D([0,1],[0,0],[0,0],linewidth=2,color=cmap(k))[0] for k in range(M.shape[1])]
        hdla=ax.plot3D([0,1],[0,0],[0,0],linewidth=3,color='black')[0]
        
        for M0 in M:
            self.single_frame()
            for M00,h in zip(M0,hdl):
                h.set_data_3d([0,M00[0]],[0,M00[1]],[0,M00[2]])
            MA=M0.mean(0)
            hdla.set_data_3d([0,MA[0]],[0,MA[1]],[0,MA[2]])
        self.single_frame()
        
        return v
        
    def movie_hahn_dynamic(self,v:list=None,ax=None):
        """
        Generates stochastic trajectories and plots these, including a pi pulse
        to refocus them

        Parameters
        ----------
        ax : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        
        v=np.random.rand(11)*(self.v2-self.v1)+self.v1 if v is None else np.array(v)
        n=len(v)
        v1,v2=self.v1,self.v2
        dv=(v1-v2)*np.array([-0.5,0.5])
        
        if ax is None:
            self.fig=plt.figure()
            ax=self.fig.add_subplot(1,1,1,projection='3d')
        
        self.make3Daxis(ax)
        nt,tf=225,1/v.mean()
        
        c1,s1=[getattr(np,x)(np.linspace(0,np.pi,10)) for x in ['cos','sin']]
        Mall=list()
        for k,v0 in enumerate(v):
            self.v1,self.v2=v0+dv
            out=self.single_spin_traj(state0=1 if np.mod(k,2) else 2,nt=nt,tf=tf)
            M=np.array([out['I'].real,out['I'].imag,np.zeros(out['I'].shape)])
            M1=np.array([np.ones(c1.shape)*M[0][-1],c1*M[1][-1],-s1*M[1][-1]])       
            out=self.single_spin_traj(state0=out['state'][-1],nt=nt,tf=tf)
            M2=np.array([out['I'].real*M1[0][-1]-out['I'].imag*M1[1][-1],
                         out['I'].imag*M1[0][-1]+out['I'].real*M1[1][-1],M1[2][-1]*np.ones(out['I'].shape)])
            Mall.append(np.concatenate((M,M1,M2),axis=1))
        Mall=np.moveaxis(Mall,-1,0)
        
        cmap=plt.get_cmap('tab10')
        hdl=[ax.plot3D([0,1],[0,0],[0,0],linewidth=2,color=cmap(k))[0] for k in range(n)]
        hdla=ax.plot3D([0,1],[0,0],[0,0],linewidth=3,color='black')[0]
        
        for M0 in Mall:
            self.single_frame()
            for M00,h in zip(M0,hdl):
                h.set_data_3d([0,M00[0]],[0,M00[1]],[0,M00[2]])
            MA=M0.mean(0)
            hdla.set_data_3d([0,MA[0]],[0,MA[1]],[0,MA[2]])
        self.single_frame()
        
        self.v1,self.v2=v1,v2
        
        
        
        
        
        
        
        

    #%% Movie management tools                
    def save_movie(self):
        """
        Final step in saving a movie

        Returns
        -------
        None.

        """
        
        frame_rate=self.fr
        
        filename=os.path.join(self._moviefilename.split('.')[0],'frame{0}.png')
        img0=cv2.imread(filename.format(0))
        height,width,layers=img0.shape
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video=cv2.VideoWriter(self._moviefilename,fourcc,frame_rate,(width,height))
        
        k=0
        while True:
            if os.path.exists(filename.format(k)):
                video.write(cv2.imread(filename.format(k)))
                os.remove(filename.format(k))
                k+=1
            else:
                break
        cv2.destroyAllWindows()
        video.release()
        os.rmdir(self._moviefilename.split('.')[0])
        self._moviefilename=None
        

    def make_movie(self,filename:str,function:str,**kwargs):
        """
        Make a movie out of the function "function" and save the result.

        Parameters
        ----------
        filename : str
            Name to save movie under.
        function : str
            Function to call for the movie.
        **kwargs : TYPE
            Arguments to pass for the movie.

        Returns
        -------
        None.

        """
        
        
        
        self.movie_reset()
        self._moviefilename=filename
        getattr(self,function)(**kwargs)
        self.save_movie()
        self._moviefilename=None
        
        

    def movie_reset(self):
        """
        Resets the frame counter for movies.

        Returns
        -------
        None.

        """
        self._frind=0
        
        
    def single_frame(self):
        """
        Either pauses to display the current frame, or saves the current frame
        as a png for building into a movie.
        

        Parameters
        ----------
        filename : str, optional
            Filename if a movie should be saved. The default is None.
        reset : bool, optional
            Run this to set the frame counter back to 0. The default is False.

        Returns
        -------
        None.

        """
        frame_rate=self.fr
        
        filename=self._moviefilename
        if filename is None:
            plt.pause(1/frame_rate)
            self.fig.canvas.draw()
        else:
            if not(os.path.exists(filename.split('.')[0])):
                os.mkdir(filename.split('.')[0])
            self.fig.savefig(os.path.join(filename.split('.')[0],'frame{0}.png'.format(self._frind)))
            self._frind+=1
    
    
        
        
        
        