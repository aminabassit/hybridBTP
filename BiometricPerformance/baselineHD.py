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
    
    
    dataset = 'IITD'
    
    dataDir = './data/IITD_IrisCodes/lg_'
    datasetDir = dataset+'_lg'
    subjects = 224
    subjectIDs = list(range(1, subjects+1))
    
    
    precisionD = 1E-6
    

    print('Hamming Distance on dataset = ', dataset)


    

    
    print('-- Run mated and non-mated comparison for Hamming Distance')
        
    pool = mp.Pool(32)
    
    
    gen_args = zip(repeat(dataDir), subjectIDs)

    imp_args = zip(repeat(dataDir), subjectIDs)
    
    print('Start pool.starmap')
    print('Hamming distance mated')
    
    matedSc = pool.starmap(matedHD_IITD, gen_args)
    gc.collect()
    
    print('Hamming distance non-mated')

    nonMatedSC = pool.starmap(nonMatedHD_IITD, imp_args)
    gc.collect()
    
    print('End pool.starmap')
    
    pool.close()
    pool.join()
    
    matedSc = flattenList(matedSc)
    nonMatedSC = flattenList(nonMatedSC)
    
    print('len(matedSc) = ', len(matedSc), ' len(nonMatedSC) = ', len(nonMatedSC))
    print('matedSc = ', matedSc[:10], ' nonMatedSC = ', nonMatedSC[:10])
    
    
    
    print('-- Measure the performance')
    
    resultsFile = './results/{}/Results_HD_{}_{}.pkl'.format(datasetDir, dataset,timeTag)
    
    
    plotFile = './results/{}/DET_HD_{}_{}'.format(datasetDir, dataset,timeTag)
    plotTitle = 'HD DET curve of {} dataset'.format(dataset)  

    

    _, fmr, tmr = get_fmr_tmr_Prec(-np.array(matedSc), -np.array(nonMatedSC), precision = precisionD)
    eer = calc_eerPT(fmr, tmr)
    print('eer = ', eer)
    gc.collect() 
    pickle.dump((matedSc, nonMatedSC, fmr, tmr, eer), open(resultsFile, 'wb'))

    fmrList = [fmr]
    tmrList = [tmr] 
    labelList = ['Hamming Distance (EER = {:.2E})'.format(eer)] 
    linewidthList = [2]
    
    plot_MultipleDETsSaved(plotFile+'.png', plotTitle, fmrList, tmrList, labelList, linewidthList)

    



if __name__ == '__main__':
    main()












