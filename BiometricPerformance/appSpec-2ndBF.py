#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from utils import *
from bloomFilter import *

import gc
import pickle
import numpy as np
from itertools import repeat
import multiprocessing as mp
from datetime import datetime






def main():    
    now = datetime.now()
    timeTag = now.strftime("%d%m%Y_%H%M%S") 

    
    hb = 10
    hB = 20
    wb = 32     
    wB = 512
    
    nBlocks = (wB//wb) * (hB//hb)
   
    
    

    dataset = 'IITD'
    dataDir = './data/IITD_IrisCodes/lg_'
    datasetDir = dataset+'_lg_AppSpec'
    matedHD = mated2ndBF_IITD
    nonMatedHD = nonMated2ndBF_IITD
    subjects = 224
    subjectIDs = list(range(1, subjects+1))
     


    
    precisionD = 1E-6
    

    print('2nd BF on dataset = ', dataset, ' Application Specific Setting')


    permKey = getPermutationKey(80)  
    permutFile = './results/{}/permKey_2ndBF_HD_{}_{}.pkl'.format(datasetDir, dataset,timeTag)
    pickle.dump(permKey, open( permutFile, 'wb'))

    
    print('-- Run mated and non-mated comparison for 2nd BF Hamming Distance')
        
    pool = mp.Pool(64)
    
    
    gen_args = zip(repeat(dataDir), subjectIDs, repeat(hb), repeat(wb), repeat(permKey))

    imp_args = zip(repeat(dataDir), subjectIDs, repeat(hb), repeat(wb), repeat(permKey))
    
    print('Start pool.starmap')
    print('Hamming distance mated')
    
    matedSc = pool.starmap(matedHD, gen_args)
    gc.collect()
    
    print('Hamming distance non-mated')

    nonMatedSC = pool.starmap(nonMatedHD, imp_args)
    gc.collect()
    
    print('End pool.starmap')
    
    pool.close()
    pool.join()
    
    matedSc = flattenList(matedSc)
    nonMatedSC = flattenList(nonMatedSC)
    
    print('matedSc = ', len(matedSc), ' nonMatedSC = ', len(nonMatedSC))
    print('matedSc = ', matedSc[:5], ' nonMatedSC = ', nonMatedSC[:5])

    
    
    print('-- Measure the performance')
    
    _, fmr, tmr = get_fmr_tmr_Prec(-np.array(matedSc), -np.array(nonMatedSC), precision = precisionD)
    eer = calc_eerPT(fmr, tmr)
    print('eer = ', eer)
    gc.collect() 
    
    resultsFile = './results/{}/Results_2ndBF_HD_{}_{}.pkl'.format(datasetDir, dataset,timeTag)
    pickle.dump((matedSc, nonMatedSC, fmr, tmr, eer), open(resultsFile, 'wb'))
    
    plotFile = './results/{}/DET_2ndBF_HD_{}_{}.png'.format(datasetDir, dataset,timeTag)
    plotTitle = '2nd BF HD DET curve of {} dataset'.format(dataset)  

    fmrList = [fmr]
    tmrList = [tmr] 
    labelList = ['Hamming Distance (EER = {:.2E})'.format(eer)] 
    linewidthList = [2]
    
    plot_MultipleDETsSaved(plotFile, plotTitle, fmrList, tmrList, labelList, linewidthList)

    










if __name__ == '__main__':
    main()


