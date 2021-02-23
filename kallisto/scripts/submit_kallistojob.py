#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 23:31:15 2019

@author: sturkars
"""
import glob
import sys
import os
import re

# data and results directories
org = "R12"
run_dir = "/proj/omics4tb2/ENIGMA/Baliga_pH"
data_dir = "%s/RawData/20201223_Discovery_ENIGMA" %run_dir
data_folders = glob.glob('%s/*' %(data_dir))
data_folders = [element for element in data_folders if element not in ('%s/*/trimmed,%s/*/fastqc_results,%s/EPD3*')%(data_dir,data_dir,data_dir)]

#print(data_folders)
print('Total Samples to process: %s' %(len(data_folders)))
#sys.exit()

jobscripts_dir = "%s/scripts/kallisto_jobscripts_%s" %(run_dir,org)
jobscripts_logs = "%s/kallisto_logs_%s" %(jobscripts_dir,org)

# create sample specific results directory
if not os.path.exists('%s' %(jobscripts_dir)):
    os.makedirs('%s' %(jobscripts_dir))

if not os.path.exists('%s' %(jobscripts_logs)):
    os.makedirs('%s' %(jobscripts_logs))

folderCount = 1
for data_folder in data_folders:
    folder_name = data_folder.split('/')[-1]
    job_name = 'j-%s' %folder_name
    jobscript = '%s/%s_%s_kallisto.csh' %(jobscripts_dir,folder_name,org)
    # write to job file
    with open(jobscript,'w') as g:
      g.write('#!/bin/bash\n\n')
      g.write('#$ -N %s\n'%(job_name))
      g.write('#$ -o %s/%s_outlog.txt\n' %(jobscripts_logs,job_name))
      g.write('#$ -e %s/%s_errorlog.txt\n' %(jobscripts_logs,job_name))
      g.write('#$ -pe smp 2\n')
      g.write('#$ -S /bin/bash\n\n')

      g.write('#Sample: %s\n' %(data_folder))

      # changing terminal to bash
      g.write('bash\n\n')
      g.write('source /users/sturkars/.bashrc \n\n')

      # change directory
      g.write('cd %s/scripts\n\n' %(run_dir))

      job_cmd = 'python run_kallisto_%s.py %s' %(org,folder_name)
      print(job_cmd)
      # spades command
      g.write('%s' %(job_cmd))

    g.close()
    folderCount = folderCount + 1

    # submit each job with qsub
    cmd = 'qsub %s' %jobscript
    print(cmd)
    print
    os.system(cmd)
    #sys.exit()

