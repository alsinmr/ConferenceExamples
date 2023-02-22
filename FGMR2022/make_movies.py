#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 10:19:37 2022

@author: albertsmith
"""

from Exchange import Exchange
from T1movie import T1movie
import numpy as np
import matplotlib.pyplot as plt

#%% Create Exchange movies

ex=Exchange()
ex.fr=30
ex.make_movie('single_spin1.mp4','single_spin_v_rate',state0=2)
ex.fr=15
ex.make_movie('accum.mp4')

fig=plt.figure()
ax=[fig.add_subplot(1,3,k+1) for k in range(3)]
titles=[r'Slow ($k=1/10*k_{coalescence}$)',r'Intermediate ($k=k_{coalescence}$)',r'Fast ($k=10*k_{coalescence}$)']
# ex.nt=
for k,a,title in zip(ex.kcoal*np.array([1/ex.rate_scaling,1,ex.rate_scaling]),ax,titles):
    ex.k12=k
    out=ex.calc_traj(tf=5)
    a.plot(out['f'],out['S'])
    a.set_xlabel(r'$\delta$ / Hz')
    a.set_title(title)
    a.set_xlim([-20,30])
    a.set_ylim([0,250])
fig.set_size_inches([9.4,4])
fig.tight_layout()

ex.make_movie('ksweep.mp4','movie_ksweep',ratio=np.logspace(-1.5,1.5,150),tf=3)

ex.fr=30
v=np.random.rand(10)*5+5
ex.make_movie('hahn_echo_static.mp4','movie_hahn_static',v=v)
ex.v1=5
ex.v2=10
ex.k12=50
ex.make_movie('hahn_echo_dynamic.mp4','movie_hahn_dynamic',v=v)

#%% Create T1 movies
T1=T1movie(v0=2,CSA=.4,delta=0.1)
T1.nt=900

T1.nsteps=500
T1.make_movie('T1fast.mp4','movie_all')

while T1.M[:,2].max()<.75:T1.nsteps=14
T1.make_movie('T1match.mp4','movie_all')

while T1.M[:,2].max()>-.9:T1.nsteps=1
T1.make_movie('T1slow.mp4','movie_all')

T1.nsteps=14
T1.make_movie('T1accum.mp4','movie_accu')