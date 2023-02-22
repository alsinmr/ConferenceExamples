#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 12:18:30 2022

@author: albertsmith
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import cv2

class T1movie():
    def __init__(self,v0=2,theta=np.pi,CSA=.4,nsteps=14,delta=0.1,dt=1/30,nt=450,fr=30):
        #Fast: nsteps=500; Match: nsteps=14; slow: nsteps=1
        self.v0=v0
        self.theta=theta
        self.delta=delta
        self.nsteps=nsteps
        self.dt=dt
        self.nt=nt
        self.fr=fr
        self.CSA=CSA
        self._vCSA=None
        self._M=None
        self._moviefilename=None
        self.fig=None
        self.relax2eq=True
    
    def __setattr__(self, name, value):
        if name in ['delta','nt','nsteps']:
            self._vCSA=None
        if name in ['delta','nt','theta','CSA','v0','dt','_vCSA']:
            self._M=None
        super().__setattr__(name, value)
    
    @property
    def t(self):
        """
        Return time axis

        Returns
        -------
        np.array
            Time axis.

        """
        return np.arange(self.nt*self.nsteps)*self.dt/self.nsteps
    
    @property
    def vCSA(self):
        """
        Returns the direction of the CSA

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if self._vCSA is None:
            return self.gen_traj()
        return self._vCSA
    

    @property
    def Ct(self):
        """
        Calculates the correlation function corresponding to the trajectory of
        the CSA

        Returns
        -------
        None.

        """
        Ct=np.zeros(self.vCSA.shape[0]+1,dtype=complex)
        for k in range(3):
            for j in range(k,3):
                FT=np.fft.rfft(self.vCSA[:,k]*self.vCSA[:,j],n=self.vCSA.shape[0]<<1)
                Ct+=FT.conj()*FT*(3/2 if k==j else 3)
        Ct=np.fft.irfft(Ct)[:self.vCSA.shape[0]]
        Ct/=np.arange(len(Ct),0,-1)
        Ct-=1/2
        return Ct
    
    @property
    def S2(self):
        """
        Estimates S2 for the trajectory of the CSA

        Returns
        -------
        None.

        """
        
        S2=-1/2
        for k in range(3):
            for j in range(k,3):
                S2+=(self.vCSA[:,k]*self.vCSA[:,j]).mean()**2*(3/2 if k==j else 3)
        return S2
    
    @property
    def tc(self):
        """
        Estimates the correlation time of the motion

        Returns
        -------
        None.

        """
        Ct=(self.Ct-self.S2)/(1-self.S2)
        i=np.argwhere(Ct<1-np.exp(-1))[0,0]
        return self.t[i]
                
    
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
        
    
    def gen_traj(self):
        """
        Generates and stores a trajectory for the CSA tensor (stored as a vector).
        Re-run to generate a new trajectory.

        Returns
        -------
        np.array
            Nx3 trajectory

        """
        
        v=np.zeros([self.nt*self.nsteps,3])
        v[0]=[0,0,1]
        delta=self.delta
        # delta=0.1
        
        vc=v[0]
        for k in range(self.nt*self.nsteps-1):
            for _ in range(1):
                ct=np.min([vc[2],1])
                st=np.sqrt(1-ct**2)
                cp,sp=vc[:2]/np.sqrt((vc[:2]**2).sum()) if ct<1 else (1,0)
                    
                ra=np.random.rand(1)[0]*2*np.pi
                for _ in range(2):                
                    "Frame of v[k]–> vector tilted from z with angle delta in a random (ra) direction"
                    vn0=np.array([np.sin(delta)*np.cos(ra),np.sin(delta)*np.sin(ra),np.cos(delta)])
                    
                    vn=np.array([ct*vn0[0]+st*vn0[2],vn0[1],-st*vn0[0]+ct*vn0[2]])
                    vn=np.array([cp*vn[0]-sp*vn[1],sp*vn[0]+cp*vn[1],vn[2]])
                    
                    if vn[2]>=np.cos(self.theta):break
                    ra+=np.pi  #Flip the angle 180 degrees if we go over the opening cone angle
                vc=vn #Update vc
            v[k+1]=vn
            
        self._vCSA=v
        return v
    
    @property
    def M(self):
        """
        Returns the direction of the magnetization as a function of time

        Returns
        -------
        None.

        """
        
        if self._M is None:
            M=np.zeros([self.nt,3])
            M[0]=[0,0,-1]
            Mc=M[0]
            counter=1
            for k,vCSA in enumerate(self.vCSA[:-1]):
                
                ct=np.min([vCSA[2],1])
                st=np.sqrt(1-ct**2)
                s2t=2*ct*st
                cp,sp=vCSA[:2]/np.sqrt((vCSA[:2]**2).sum()) if ct<1 else (1,0)
                
                
                
                v=self.CSA*np.array([-3/4*s2t*cp,3/4*s2t*sp,(3*ct**2-1)/2])
                v[2]+=self.v0
                
                # v=self.CSA*vCSA+np.array([0,0,self.v0])  #CSA plus Larmor frequency
                
                veff=np.sqrt((v**2).sum())
                vn=v/veff #Normalized v
                veff*=2*np.pi
                
                ct=np.min([vn[2],1])
                st=np.sqrt(1-ct**2)
                cp,sp=vn[:2]/np.sqrt((vn[:2]**2).sum()) if ct<1 else (1,0)
                
                for _ in range(1+self.relax2eq):
                    "Bring into frame of interactions"
                    M0=np.array([cp*Mc[0]+sp*Mc[1],-sp*Mc[0]+cp*Mc[1],Mc[2]])
                    M0=np.array([ct*M0[0]-st*M0[2],M0[1],st*M0[0]+ct*M0[2]])
                    "Rotate about interaction"
                    ci,si=np.cos(veff*self.dt/self.nsteps),np.sin(veff*self.dt/self.nsteps)
                    M0=np.array([ci*M0[0]-si*M0[1],si*M0[0]+ci*M0[1],M0[2]])
                    "Rotate back into lab frame"
                    M0=np.array([ct*M0[0]+st*M0[2],M0[1],-st*M0[0]+ct*M0[2]])
                    M0=np.array([cp*M0[0]-sp*M0[1],sp*M0[0]+cp*M0[1],M0[2]])
                    # veff*=1+0.25*(Mc[2]<M0[2])
                
                Mc=M0
                if np.mod(k+1,self.nsteps)==0:
                    M[counter]=M0
                    counter+=1
            self._M=M
                
        return self._M
    
    #%% Movies
    def movie_CSA(self,ax=None):
        if ax is None:
            ax=plt.figure().add_subplot(111,projection='3d')
        self.make3Daxis(ax)
        
        hdl=ax.plot3D([0,self.vCSA[0][0]],[0,self.vCSA[0][1]],[0,self.vCSA[0][2]],color='magenta',linewidth=3)[0]
        # return hdl
        for v in self.vCSA[::self.nsteps]:
            self.single_frame()
            hdl.set_data_3d([0,v[0]],[0,v[1]],[0,v[2]])
        
    def movie_M(self,ax=None):
        if ax is None:
            ax=plt.figure().add_subplot(111,projection='3d')
        self.make3Daxis(ax)
        
        hdl=ax.plot3D([0,self.M[0][0]],[0,self.M[0][1]],[0,self.M[0][2]],color='black',linewidth=3)[0]
        # return hdl
        for v in self.M:
            self.single_frame()
            hdl.set_data_3d([0,v[0]],[0,v[1]],[0,v[2]])
            
    def movie_all(self):
        self.fig=plt.figure()
        self.fig.set_size_inches([10,7.25])
        ax=[self.fig.add_subplot(2,2,k+1,projection='3d') for k in range(4)]
        
        titles=['CSA (lab)','CSA (rotating)','<M> (lab)','<M> (rotating)']
        hdl=list()
        for k,(a,title) in enumerate(zip(ax,titles)):
            self.make3Daxis(a)
            a.set_title(title)
            hdl.append(a.plot3D([0,0],[0,0],[0,1],color='magenta' if k<2 else 'black',linewidth=3)[0])
            
        for a in ax[:2]:
            a.view_init(azim=45,elev=90)
        self.fig.tight_layout()
        
        for k,(v,M) in enumerate(zip(self.vCSA[::self.nsteps],self.M)):
            self.single_frame()
            c,s=np.cos(k*self.dt*self.v0*2*np.pi),np.sin(k*self.dt*self.v0*2*np.pi)
            vr=[v[0]*c+v[1]*s,-v[0]*s+v[1]*c,v[2]]
            Mr=[M[0]*c+M[1]*s,-M[0]*s+M[1]*c,M[2]]
            hdl[0].set_data_3d([0,v[0]],[0,v[1]],[0,v[2]])
            hdl[1].set_data_3d([0,vr[0]],[0,vr[1]],[0,vr[2]])
            hdl[2].set_data_3d([0,M[0]],[0,M[1]],[0,M[2]])
            hdl[3].set_data_3d([0,Mr[0]],[0,Mr[1]],[0,Mr[2]])

    def movie_accu(self,nspins=20):
        
        nsteps=self.nsteps
        
        M0=list()
        for k in range(nspins):
            M0.append(self.M)
            self.nsteps=nsteps #This resets the trajectory
        
        self.fig=plt.figure()
        self.fig.set_size_inches([9.5,4.8])
        ax=[self.fig.add_subplot(1,2,k+1,projection='3d') for k in range(2)]
        
        titles=['<M> (lab)','<M> (rotating)']
        hdl=list()
        hdlr=list()
        for k,(a,title,hdl0) in enumerate(zip(ax,titles,(hdl,hdlr))):
            self.make3Daxis(a)
            a.set_title(title)
            for _ in range(nspins):
                hdl0.append(a.plot3D([0,0],[0,0],[0,1],color='cyan',alpha=.25,linewidth=2)[0])
        hdla=ax[0].plot3D([0,0],[0,0],[0,1],color='black',linewidth=3)[0]
        hdlar=ax[1].plot3D([0,0],[0,0],[0,1],color='black',linewidth=3)[0]
        
        for k,M in enumerate(zip(*M0)): #This sweeps over time, yielding the direction of each magnetization vector
            self.single_frame()
            c,s=np.cos(k*self.dt*self.v0*2*np.pi),np.sin(k*self.dt*self.v0*2*np.pi)
            for h,hr,M1 in zip(hdl,hdlr,M):
                h.set_data_3d([0,M1[0]],[0,M1[1]],[0,M1[2]])
                Mr=[M1[0]*c+M1[1]*s,-M1[0]*s+M1[1]*c,M1[2]]
                hr.set_data_3d([0,Mr[0]],[0,Mr[1]],[0,Mr[2]])
            Ma=np.mean(M,axis=0)
            hdla.set_data_3d([0,Ma[0]],[0,Ma[1]],[0,Ma[2]])
            Mra=[Ma[0]*c+Ma[1]*s,-Ma[0]*s+Ma[1]*c,Ma[2]]
            hdlar.set_data_3d([0,Mra[0]],[0,Mra[1]],[0,Mra[2]])
        
            
        
        
        
        
    
    #%% Movie management tools                
    def save_movie(self):
        """
        Final step in saving a movie

        Returns
        -------
        None.

        """
        filename=os.path.join(self._moviefilename.split('.')[0],'frame{0}.png')
        img0=cv2.imread(filename.format(0))
        height,width,layers=img0.shape
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video=cv2.VideoWriter(self._moviefilename,fourcc,self.fr,(width,height))
        
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
        if filename is not None:plt.close(self.fig) #Faster?
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
        filename=self._moviefilename
        if filename is None:
            plt.pause(1/self.fr)
            self.fig.canvas.draw()
        else:
            if not(os.path.exists(filename.split('.')[0])):
                os.mkdir(filename.split('.')[0])
            self.fig.savefig(os.path.join(filename.split('.')[0],'frame{0}.png'.format(self._frind)))
            self._frind+=1
    
        
        
        
        
    