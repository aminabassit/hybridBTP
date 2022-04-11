#!/usr/bin/env python3
# -*- coding: utf-8 -*-



from utils import *
from bloomFilter import *

import gc
import pickle
import numpy as np

import glob 
import os

    


hb = 10
hB = 20
wb = 32     
wB = 512

nBlocks = (wB//wb) * (hB//hb)
   



dataset = 'IITD'



dataDir = './data/IITD_IrisCodes/lg_' 



dataDest = './data/ITD_IrisCodes_left_lines/'

for sub in range(1,225):
    for emb in range(1,6):
        nameF = dataDir+int2str(sub,3)+int2str(emb,2)+'.png.txt'
        x = readIrisCode(nameF)
        xB = x
        nameFDest = 'lg_'+str(sub)+'_'+str(emb)+'.txt'
        np.savetxt(dataDest+nameFDest, xB, delimiter=',',fmt='%d' )


dataDest = './data/IITD_IrisCodes_left_blocks/' 
for sub in range(1,225):
    for emb in range(1,6):
        nameF = dataDir+int2str(sub,3)+int2str(emb,2)+'.png.txt'
        x = readIrisCode(nameF)
        xB = splitIrisCodeToBlocks(hb, wb, x)
        xBT = [xB[i].T for i in range(len(xB))]
        xBT = np.concatenate(xBT)
        nameFDest = 'lg_'+str(sub)+'_'+str(emb)+'.txt'
        np.savetxt(dataDest+nameFDest, xBT, delimiter=',',fmt='%d' )




















