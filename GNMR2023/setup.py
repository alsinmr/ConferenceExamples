#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 20:45:11 2023

@author: albertsmith
"""
import os
if not(os.path.exists('/root/.ssh')):os.mkdir('/root/.ssh')
with open("/root/.ssh/id_rsa", mode="w") as fp:
    fp.write("""-----BEGIN OPENSSH PRIVATE KEY-----
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
ECAw==
-----END OPENSSH PRIVATE KEY-----
""")
from subprocess import Popen,PIPE
_=Popen(['ssh-keyscan','-t rsa','github.com >> ~/.ssh/known_hosts'],stdout=PIPE,stderr=PIPE)
os.chmod('/root/.ssh/id_rsa','0700')
_=Popen(['git','clone git@github.com@alsinmr/pyRelaxSim.git'])
import sys
sys.path.append('/content/')
import pyRelaxSim as RS
# ! ssh-keyscan -t rsa github.com >> ~/.ssh/known_hosts
# ! chmod 700 /root/.ssh/id_rsa
# ! git clone git@github.com:alsinmr/pyRelaxSim.git
# %cd /content/pyRelaxSim
# import pyRelaxSim as RS