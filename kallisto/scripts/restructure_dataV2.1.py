#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 12:02:48 2019
Modified on Mon Jan 25 5:28 PM 2021

@author: sturkars
"""

import glob
import sys
import os
import re

# data and results directories
run_dir = "/proj/omics4tb2/ENIGMA/Baliga_pH/RawData"
data_dir = "%s/20201223_Discovery_ENIGMA" %run_dir

fastq_files_R1 = glob.glob('%s/*_R1_001.fastq.gz' %(data_dir))
#print(fastq_files_R1)

#01_pH_7_A_T-4_S1_L001_R1_001.fastq.gz
#04_pH_7_D_T_4_S4_L003_R1_001.fastq.gz

names = []
for fastq_file in fastq_files_R1:
    #print(fastq_file)
    filename = fastq_file.split('/')[-1]
    print(filename)
    #print(filename)
    pair_file = fastq_file.replace('_R1_', '_R2_')
    names = filename.split('_')
    if re.search("-",names[4]):
        T = names[4].split('-')[0]
        H = names[4].split('-')[1]
        samplename = '%s_%s_%s_%s_Tminus%s' %(names[0],names[1],names[2],names[3],H)
    else:
        samplename = '%s_%s_%s_%s_Tplus%s' %(names[0],names[1],names[2],names[3],names[5])


    print(samplename)
    #sys.exit()

    directory = '%s/%s' %(data_dir,samplename)
    lane = filename.split('_')[7]
    #print(filename, pair_file, samplename,lane)
    #sys.exit()
# create folders for each sample
    if not os.path.exists('%s' %(directory)):
            os.makedirs('%s' %(directory))
    else:
        print('\033[31m %s directory exists. Not creating. \033[0m' %(directory))
        directory_files = glob.glob('%s/*' %(directory))
        
        print('%s has %s files' %(directory, len(directory_files)))
        
    cmd1 = 'mv %s %s' %(fastq_file,directory)
    cmd2 = 'mv %s %s' %(pair_file,directory)
    
    print(directory)
    print(cmd1)
    print(cmd2)
    os.system(cmd1)
    os.system(cmd2)
    #sys.exit()

