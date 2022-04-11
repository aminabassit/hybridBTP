#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 22:02:47 2022

@author: Amina
"""




import numpy as np
import glob
import pickle

from utils import *
from scipy.spatial import distance


"""
Functions in common

"""

def genBloomFilter(transfBlock):
    n,m = transfBlock.shape
    bfTemp = np.zeros(2**n, dtype=int)    
    indOnes = list(set(binArray2int(transfBlock[:,i]) for i in range(m)))
    bfTemp[indOnes] = [1]*len(indOnes)    
    return bfTemp


def genBFbasedTemplate(transBlocks):
    BFs = [genBloomFilter(transB) for transB in transBlocks]
    return np.vstack(BFs).T






"""
Functions specific to BF 2nd category

"""




# n = (hb * wB // (sqrt(regions) * wb) ) 
def getPermutationKey(n, regions = 4):
    permutation = [ np.random.permutation(n) for i in range(regions)]
    return np.asarray(permutation)

def permuteWithinRegion(hb, wb, features, permutation, regions = 4):
    m,n = features.shape
    r = int(np.sqrt(regions))
    featRegions = [f.reshape((hb * n//(wb*r), wb)) for feat in np.hsplit(features, r) for f in np.vsplit(feat, r)]    
    featPerm = [f[permutation[i]].reshape((hb, wb * n//(wb*r) )) for i, f in enumerate(featRegions)] 
    featPermBlocks = [b for block in featPerm for b in np.hsplit(block, block.shape[1]//wb)]   
    return featPermBlocks




def matedHD_IITD(dataDir, subject): 
    
    matedSc = []
    
    subjectEmbs = glob.glob(dataDir+int2str(subject,3)+'*.png.txt')
    subjectEmbs.sort()
    
    
    
    embs = [ emb for (i, emb) in enumerate(subjectEmbs) if i in range(0,5)]
    
    
    indexes = [(i,j) for i in range(len(embs)) for j in range(len(embs)) if i<j]
    
    
    for ind in indexes: 
        
        (tempInd, probeInd) = ind
        
        temp = readIrisCode(embs[tempInd])
        probe = readIrisCode(embs[probeInd])     
        matedSc.append(hammingIrisCodeShift(temp, probe, 8))
        
    
    return matedSc



def nonMatedHD_IITD(dataDir, subject): 
    
    nonMatedSc = []
    subjectList = list(range(1, 225))
    subjectList = subjectList[subject:]
    
    lefts = list(range(1,6))    
    fixIR = [dataDir+int2str(subject,3)+int2str(i,2)+'.png.txt' for i in lefts]
    fixIR.sort()
    impIR = [ dataDir+int2str(subj,3)+int2str(2,2)+'.png.txt' for subj in subjectList]   
    
    for tempIR in fixIR:
        temp = readIrisCode(tempIR)  
        for probeIR in impIR: 
            probe = readIrisCode(probeIR)     
            nonMatedSc.append(hammingIrisCodeShift(temp, probe, 8))               
   
         
    return nonMatedSc




def mated2ndBF_IITD(dataDir, subject, hb, wb, permutation): 
    
    matedSc = []
    
    subjectImgs = glob.glob(dataDir+int2str(subject,3)+'*.png.txt')
    subjectImgs.sort()
    subjectImgs = subjectImgs[:5]
    
    
    indexes = [(i,j) for i in range(len(subjectImgs)) for j in range(len(subjectImgs)) if i<j]
    for ind in indexes:        
        (i,j) = ind
        x = readIrisCode(subjectImgs[i])
        xBlocks = permuteWithinRegion(hb, wb, x, permutation)
        xBF = genBFbasedTemplate(xBlocks)    
        
        y = readIrisCode(subjectImgs[j])
        yBlocks = permuteWithinRegion(hb, wb, y, permutation)
        yBF = genBFbasedTemplate(yBlocks)
        
        matedSc.append(whDistBFTemplate(xBF, yBF))
    
    return matedSc



def nonMated2ndBF_IITD(dataDir, subject, hb, wb, permutation): 
    
    nonMatedSc = []  
    
    subjectList = list(range(1, 225))
    subjectList = subjectList[subject:]
    
    lefts = list(range(1,6))    
    fixIR = [dataDir+int2str(subject,3)+int2str(i,2)+'.png.txt' for i in lefts]
    fixIR.sort()
    impIR = [ dataDir+int2str(subj,3)+int2str(2,2)+'.png.txt' for subj in subjectList]
    
    for tempIR in fixIR:
        x = readIrisCode(tempIR)   
        xBlocks = permuteWithinRegion(hb, wb, x, permutation)
        xBF = genBFbasedTemplate(xBlocks) 
        for probeIR in impIR: 
            y = readIrisCode(probeIR) 
            yBlocks = permuteWithinRegion(hb, wb, y, permutation)
            yBF = genBFbasedTemplate(yBlocks)            
            nonMatedSc.append(whDistBFTemplate(xBF, yBF))
    
    return nonMatedSc


def mated2ndBF_IITD_diff(dataDir, subject, hb, wb, permutation): 
    
    matedSc = []
    
    subjectImgs = glob.glob(dataDir+int2str(subject,3)+'*.png.txt')
    subjectImgs.sort()
    subjectImgs = subjectImgs[:5]
    
    
    indexes = [(i,j) for i in range(len(subjectImgs)) for j in range(len(subjectImgs)) if i<j]
    for ind in indexes:        
        (i,j) = ind
        x = readIrisCode(subjectImgs[i])
        xBlocks = permuteWithinRegion(hb, wb, x, permutation)
        xBF = genBFbasedTemplate(xBlocks)    
        
        y = readIrisCode(subjectImgs[j])
        yBlocks = permuteWithinRegion(hb, wb, y, permutation)
        yBF = genBFbasedTemplate(yBlocks)
        
        matedSc.append(whDistBFTemplate(xBF, yBF))
    
    return matedSc



def nonMated2ndBF_IITD_diff(dataDir, subject, hb, wb, permutation, permutations): 
    
    nonMatedSc = []  
    
    subjectList = list(range(1, 225))
    subjectList = subjectList[subject:]
    permList = permutations[subject:]
    
    lefts = list(range(1,6))    
    fixIR = [dataDir+int2str(subject,3)+int2str(i,2)+'.png.txt' for i in lefts]
    fixIR.sort()
    impIR = [ dataDir+int2str(subj,3)+int2str(2,2)+'.png.txt' for subj in subjectList]
    impIR.sort()
    
    for tempIR in fixIR:
        x = readIrisCode(tempIR)   
        xBlocks = permuteWithinRegion(hb, wb, x, permutation)
        xBF = genBFbasedTemplate(xBlocks) 
        for probeIR, perm in zip(impIR, permList): 
            y = readIrisCode(probeIR) 
            yBlocks = permuteWithinRegion(hb, wb, y, perm)
            yBF = genBFbasedTemplate(yBlocks)            
            nonMatedSc.append(whDistBFTemplate(xBF, yBF))
    
    return nonMatedSc




"""
Functions specific to BF 1st category

"""

def genColumnKey(hb):
    return np.random.randint(2, size=(hb,1))

def genColumnKeys(hb, nBlocks):
    return [np.random.randint(2, size=(hb,1)) for i in range(nBlocks)]

def splitIrisCodeToBlocks(hb, wb, irisCode):
    blocks = []
    m,n = irisCode.shape
    irisCodeColmnSplit = np.hsplit(irisCode, n//wb)    
    blocks = [np.vsplit(b, m//hb) for b in irisCodeColmnSplit]
    return [b for block in blocks for b in block]

def getTransBlocksFirst(blocks, key, wb):
    keyB = np.tile(key, wb)
    transfBlocks = [b ^ keyB for b in blocks]
    return transfBlocks

def xorWithDiffColmKeys(blocks, key, wb):
    transfBlocks = [b ^ np.tile(k, wb) for b,k in zip(blocks, key)]
    return transfBlocks

def mated1stBF_IITD(dataDir, subject, hb, wb, key): 
    
    matedSc = []
    
    subjectImgs = glob.glob(dataDir+int2str(subject,3)+'*.png.txt')
    subjectImgs.sort()
    subjectImgs = subjectImgs[:5]    
    
    indexes = [(i,j) for i in range(len(subjectImgs)) for j in range(len(subjectImgs)) if i<j]
    for ind in indexes:        
        (i,j) = ind
        x = readIrisCode(subjectImgs[i])
        xB = splitIrisCodeToBlocks(hb, wb, x)
        # xBlocks = getTransBlocksFirst(xB, key, wb)
        xBlocks = xorWithDiffColmKeys(xB, key, wb)
        xBF = genBFbasedTemplate(xBlocks)
        
        y = readIrisCode(subjectImgs[j])
        yB = splitIrisCodeToBlocks(hb, wb, y)
        # yBlocks = getTransBlocksFirst(yB, key, wb)
        yBlocks = xorWithDiffColmKeys(yB, key, wb)        
        yBF = genBFbasedTemplate(yBlocks)
        
        matedSc.append(whDistBFTemplate(xBF, yBF))
    
    return matedSc



def nonMated1stBF_IITD(dataDir, subject, hb, wb, key): 
    
    nonMatedSc = []  
    
    subjectList = list(range(1, 225))
    subjectList = subjectList[subject:]
    
    lefts = list(range(1,6))    
    fixIR = [dataDir+int2str(subject,3)+int2str(i,2)+'.png.txt' for i in lefts]
    fixIR.sort()
    impIR = [ dataDir+int2str(subj,3)+int2str(2,2)+'.png.txt' for subj in subjectList]
    
    for tempIR in fixIR:
        x = readIrisCode(tempIR)   
        xB = splitIrisCodeToBlocks(hb, wb, x)
        # xBlocks = getTransBlocksFirst(xB, key, wb)
        xBlocks = xorWithDiffColmKeys(xB, key, wb)        
        xBF = genBFbasedTemplate(xBlocks) 
        for probeIR in impIR: 
            y = readIrisCode(probeIR) 
            yB = splitIrisCodeToBlocks(hb, wb, y)
            # yBlocks = getTransBlocksFirst(yB, key, wb)
            yBlocks = xorWithDiffColmKeys(yB, key, wb)
            yBF = genBFbasedTemplate(yBlocks)           
            nonMatedSc.append(whDistBFTemplate(xBF, yBF))
    
    return nonMatedSc



def mated1stBF_IITD_diff(dataDir, subject, hb, wb, userKey): 
    
    matedSc = []
    
    subjectImgs = glob.glob(dataDir+int2str(subject,3)+'*.png.txt')
    subjectImgs.sort()
    subjectImgs = subjectImgs[:5]    
    
    indexes = [(i,j) for i in range(len(subjectImgs)) for j in range(len(subjectImgs)) if i<j]
    for ind in indexes:        
        (i,j) = ind
        x = readIrisCode(subjectImgs[i])
        xB = splitIrisCodeToBlocks(hb, wb, x)
        xBlocks = xorWithDiffColmKeys(xB, userKey, wb)
        xBF = genBFbasedTemplate(xBlocks)
        
        y = readIrisCode(subjectImgs[j])
        yB = splitIrisCodeToBlocks(hb, wb, y)
        yBlocks = xorWithDiffColmKeys(yB, userKey, wb)        
        yBF = genBFbasedTemplate(yBlocks)
        
        matedSc.append(whDistBFTemplate(xBF, yBF))
    
    return matedSc



def nonMated1stBF_IITD_diff(dataDir, subject, hb, wb, userKey, userKeys): 
    
    nonMatedSc = []  
    
    subjectList = list(range(1, 225))
    subjectList = subjectList[subject:]
    userKeysList = userKeys[subject:]
    
    lefts = list(range(1,6))    
    fixIR = [dataDir+int2str(subject,3)+int2str(i,2)+'.png.txt' for i in lefts]
    fixIR.sort()
    impIR = [ dataDir+int2str(subj,3)+int2str(2,2)+'.png.txt' for subj in subjectList]
    impIR.sort()
    
    for tempIR in fixIR:
        x = readIrisCode(tempIR)   
        xB = splitIrisCodeToBlocks(hb, wb, x)
        xBlocks = xorWithDiffColmKeys(xB, userKey, wb)        
        xBF = genBFbasedTemplate(xBlocks) 
        for probeIR, key in zip(impIR, userKeysList): 
            y = readIrisCode(probeIR) 
            yB = splitIrisCodeToBlocks(hb, wb, y)            
            yBlocks = xorWithDiffColmKeys(yB, key, wb)
            yBF = genBFbasedTemplate(yBlocks)           
            nonMatedSc.append(whDistBFTemplate(xBF, yBF))
    
    return nonMatedSc























def matedHD_FRGC(dataDir, subject):     
    matedSc = []    
    subjectImgs = glob.glob(dataDir+int2str(subject,5)+'d*_bin.txt')    
    indexes = [(i,j) for i in range(len(subjectImgs)) for j in range(len(subjectImgs)) if i<j]
    for (i,j) in indexes:  
        x = readBinFacialFeat(subjectImgs[i])         
        y = readBinFacialFeat(subjectImgs[j])  
        matedSc.append(distance.hamming(x,y))    
    return matedSc



def nonMatedHD_FRGC(dataDir, subject, subjectIDs):
    subjIndex = subjectIDs.index(subject)    
    subjectList = subjectIDs[subjIndex+1:]    
    
    nonMatedSc = []    
    genTemp = glob.glob(dataDir+int2str(subject,5)+'d*_bin.txt')
    genTemp.sort()
    impProb = [ glob.glob(dataDir+int2str(subj,5)+'d*_bin.txt')[-1] for subj in subjectList]   
    
    for tempFR in genTemp:
        temp = readBinFacialFeat(tempFR)  
        for probeFR in impProb: 
            probe = readBinFacialFeat(probeFR)     
            nonMatedSc.append(distance.hamming(temp, probe))     
    return nonMatedSc












