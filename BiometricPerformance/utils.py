#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 21:21:04 2022

@author: Amina
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance




def flattenList(lists):
    flatList = []
    for l in lists:
        flatList += l    
    return flatList

def int2str(i, digits):
    return '0'*(digits-len(str(i)))+str(i)

def binStr2binArrayRow(binStr):
    return np.array([int(s) for s in binStr], dtype= int ).reshape((1, len(binStr)))

def binArray2int(binCol):
    binColList = list(np.array2string(binCol, formatter={'int':lambda x: str(x)}, separator=''))[1:-1]
    binColStr = ''.join(binColList)
    return int(binColStr, 2)

def readIrisCode(featFile):
    f = np.loadtxt(featFile, dtype=str).tolist()
    fList = [ binStr2binArrayRow(b) for b in f]
    iriscode = np.concatenate(fList)
    return iriscode

def readBinFacialFeat(featFile):
    f = np.loadtxt(featFile, dtype=str).tolist()
    fList = [int(b) for b in f]
    mat = np.array(fList)
    return mat


def hammingIrisCodeShift(irisCode1, irisCode2, nShift):
    hds = []
    iC1 = irisCode1.flatten()    
    for n in range(-nShift, nShift+1):
        iC2 = irisCode2
        if n !=0:
            iC2 = np.roll(iC2, n, axis=1)  
        iC2 = iC2.flatten()
        hd = np.sum(iC1 ^ iC2) / len(iC1)
        hds.append(hd)
    return min(hds)


def whDistBFTemplate(bf1, bf2):
    dist = np.sum(bf1 ^ bf2, axis = 0)
    onesBF1 = np.sum(bf1, axis = 0)
    onesBF2 = np.sum(bf2, axis = 0)
    whd = dist / (onesBF1 + onesBF2)
    return np.sum(whd)/len(dist)


# plot





def get_fmr_tmr_Prec(matedScores, nonMatedScores, precision = 0.1):
    """Calculates the FMR and TMR from mated and nonMated scores
            Parameters
            ----------
            matedScores :   array-like of shape = (n_genuine,)
                            mated scores resulting from genuine comparison of two same-subject-samples
            nonMatedScores : array-like of shape = (n_impostor,)  
                            nonMated scores resulting from impostor comparison of two different-subject-samples
            precision :     step between two adjacent threshold scores
                            equals to 1 when the scores are integers and 1E-n when they are real-values 
            
            
            Returns
            -------
            fmr :   array-like representing the False Match Rate
                    the rate of impostor scores as a function of the threshold
            tmr :   array-like representing the True Match Rate
                    the rate of genuine scores as a function of the threshold
            scores : array-like representing all possible thresholds
    """
    # determine all possible thresholds
    minS = min(min(matedScores), min(nonMatedScores))
    maxS = max(max(matedScores), max(nonMatedScores))
    scores = np.arange(minS, maxS + precision, precision)
    lenMated = len(matedScores)
    lennonMated = len(nonMatedScores)
    # calculate the occurrence of mated scores 
    pgen = np.histogram(matedScores, bins=scores, density=False)[0]
    # calculate the genuine probability 
    pgen = pgen/lenMated
    # calculate the TMR as a function of the thresholds (that is scores)
    tmr = np.maximum(0, 1 - np.cumsum(pgen))
    # calculate the occurrence of nonMated scores 
    pimp = np.histogram(nonMatedScores, bins=scores, density=False)[0]
    # calculate the impostor probability
    pimp = pimp/lennonMated
    # calculate the FMR as a function of the thresholds (that is scores)    
    fmr = np.maximum(0, 1 - np.cumsum(pimp))  
    # for a threshold equals to scores[i], FMR equals to fmr[i] and TMR equals to tmr[i] 
    return scores, fmr, tmr


def calc_eerPT(fmr, tmr):
    """Calculates the Equal Error Rate point from fmr and 1-tmr
            Parameters
            ----------
            fmr :   array-like representing the False Match Rate
                    the rate of impostor scores as a function of the threshold
            tmr :   array-like representing the True Match Rate
                    the rate of genuine scores as a function of the threshold
            
            Returns
            -------
            eer : point where fmr and fnmr (1-tmr) are equal
    """
    fnmr = 1-tmr
    x = fnmr - fmr
    if ((fmr.size == 0) or (fnmr.size == 0)):
        return np.inf
    
    index = np.argmin(np.abs(x))   
    
    if (index == 0):
        return (fnmr[index] + fmr[index])/2 

    if (index == len(fmr)):
        return (fnmr[index] + fmr[index])/2 

    if (fmr[index] <= fnmr[index]):
        l_index = index - 1
        r_index = index
    else:
        l_index = index
        r_index = index + 1

    d_1 = fmr[l_index] - fnmr[l_index]
    d_2 = fmr[r_index] - fnmr[r_index]

    if (d_1 - d_2) == 0:
        s_val = 0.5
    else:
        s_val = d_1 / (d_1 - d_2)

    eer = fnmr[l_index] + s_val*(fnmr[r_index] - fnmr[l_index])
    return eer


def plot_MultipleDETsSaved(plotFile, plotTitle, fmrList, tmrList, labelList, linewidthList):
    """Plots and saves the DET curves of fmrList and 1-tmrList
            Parameters
            ----------
            plotFile : path to where the plot should be saved
            plotTitle : path to  
            fmrList : list of FMRs 
            tmrList : list of TMRs  
            labelList : list of lables corresponding to the DET curve of the i-th fmrList and 1-tmrList
            linewidthList : list of integers specifying the linewidth of the i-th DET curve
            
            Returns
            -------
            
    """
    eer_line = np.logspace(-4,0,100) 
    fig, ax = plt.subplots()
    for fmr, tmr, lab, lw in zip(fmrList, tmrList, labelList, linewidthList):
        # fnmr is equal to 1-tmr
        ax.loglog(fmr, 1-tmr, label = lab, linewidth=lw)
    ax.loglog(eer_line, eer_line, label = 'EER')
    ax.set_aspect('equal')
    ax.set_xlim([1E-4, 1])
    ax.set_ylim([1E-4, 1])
    ax.set(xlabel='FMR', ylabel='FNMR')
    ax.legend()
    fig.suptitle(plotTitle)
    plt.savefig(plotFile)






