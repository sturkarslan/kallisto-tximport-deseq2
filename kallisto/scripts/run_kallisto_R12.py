#############################################################
##### RNASeq Analysis Pipeline with Kallisto            #####
##### Last update: 01/25/2021 Serdar Turkarslan         #####
##### Institute for Systems Biology                     #####
############################################################

import glob
import sys
import os
import re
import subprocess

# data and results directories
pipeline_version = "v.01.25.2021"
org = "R12"
run_dir = "/proj/omics4tb2/ENIGMA/Baliga_pH"
data_dir = "%s/RawData/20201223_Discovery_ENIGMA" %run_dir
genome_dir = "%s/rnaseq_reference" %run_dir
if org == "3H11":
    transcriptome_file = "%s/3H11-RAST-CDS.fa" %genome_dir
if org == "R12":
        transcriptome_file = "%s/R12-RAST-CDS.fa" %genome_dir


folder_name = str(sys.argv[1])
print(folder_name)
data_folder = "%s/%s" %(data_dir,folder_name)
#sys.exit()
############# Functions ##############

####################### Collect data files ###############################
def get_data():
    data_folders = glob.glob('%s/*' %(data_dir))
    data_folders = [element for element in data_folders if element not in ('%s/*/trimmed,%s/*/fastqc_results')%(data_dir,data_dir)]
    print('data_folders: %s, %s' %(len(data_folders),data_folders))
    return data_folders


####################### Create results directories ###############################
def create_dirs(data_trimmed_dir,fastqc_dir,results_dir): 
    dirs = [data_trimmed_dir,fastqc_dir,results_dir]
    for dir in dirs:
        # create results folder
        #print(dir)
        if not os.path.exists('%s' %(dir)):
            os.makedirs('%s' %(dir))
        else:
            print('\033[31m %s directory exists. Not creating. \033[0m' %(dir))


####################### Trimgalore for quality and trimming ###############################
def trimgalore(first_pair_file,second_pair_file,folder_name,sample_id,data_trimmed_dir,fastqc_dir):
    #print("1stpair:%s, 2ndpair:%s, folder_name:%s, sample_name:%s")%(first_pair_file,second_pair_file,folder_name,sample_name)
    print
    print ("\033[34m Running TrimGalore \033[0m")
    # create sample spepcific trimmed directory
    if not os.path.exists('%s' %(data_trimmed_dir)):
        os.makedirs('%s' %(data_trimmed_dir))
    # create sample spepcific fastqcdirectory
    if not os.path.exists('%s' %(fastqc_dir)):
        os.makedirs('%s' %(fastqc_dir))
    # run Command
    cmd = 'trim_galore --fastqc_args "--outdir %s/" --paired --output_dir %s/ %s %s' %(fastqc_dir,data_trimmed_dir,first_pair_file, second_pair_file)
    print
    print( '++++++ Trimgalore Command:', cmd)
    print
    #os.system(cmd)


####################### Collect trimmed data files ###############################
def collect_trimmed_data(data_trimmed_dir):
    # define result files
    # if file_ext == "gz":
    #     first_pair_trimmed = glob.glob('%s/*_val_1.fq.gz'%(data_trimmed_dir))
    #     second_pair_trimmed = glob.glob('%s/*_val_2.fq.gz'%(data_trimmed_dir))
    # else:
    first_pair_trimmed = glob.glob('%s/*_val_1.fq.*'%(data_trimmed_dir))
    second_pair_trimmed = glob.glob('%s/*_val_2.fq.*'%(data_trimmed_dir))
    print('Trimmed Files:\n 1st:%s \n 2nd:%s' %(first_pair_trimmed,second_pair_trimmed))
    print
    first_pair_group = ' '.join(first_pair_trimmed)
    second_pair_group = ' '.join(second_pair_trimmed)
    pair_files = []
    for file in first_pair_trimmed:
        mate_file = file.replace('R1_001_val_1.fq','R2_001_val_2.fq')
        paired_mates = file + ' ' + mate_file
        pair_files.append(paired_mates)
    
    kallisto_input_files = ' '.join(pair_files)
    
    return first_pair_group,second_pair_group, kallisto_input_files


####################### Run Salmon ###############################
def run_salmon(first_pair_group,second_pair_group,results_dir):
    print
    print('\033[33mRunning salmon! \033[0m')
    salmon_cmd = '/users/sturkars/salmon/bin/salmon quant -i %s/thaps_transcripts_index -l A -1 %s -2 %s -o %s --validateMappings --seqBias --numBootstraps 100 --gcBias' %(genome_dir, first_pair_group, second_pair_group, results_dir)
    print('Salmon run command:%s' %salmon_cmd)
    #os.system(salmon_cmd)


####################### Create Salmon index ###############################
def salmon_index():
    print
    print('\033[33mRunning salmon index! \033[0m')
    index_cmd = '/users/sturkars/salmon/bin/salmon index -t %s -i %s/thaps_transcripts_index --type quasi -k 31' %(transcriptome_file,genome_dir)
    print('salmon index command:%s' %(index_cmd))
    #os.system(index_cmd)
 
####################### Run Kalisto ###############################
def run_kallisto(first_pair_group,second_pair_group,results_dir,kallisto_input_files,org,kallisto_stdout):
    print
    print('\033[33mRunning kallisto! \033[0m')
    
    #organism specific kalisto run commands
    if org == "3H11":
        kallisto_cmd = 'kallisto quant -i %s/3H11_transcriptome_kallistoindex %s -o %s  -b 100 --bias --fusion -t 2 --rf-stranded >> %s 2>&1' %(genome_dir, kallisto_input_files, results_dir,kallisto_stdout)

    if org == "R12":
        kallisto_cmd = 'kallisto quant -i %s/R12_transcriptome_kallistoindex %s -o %s  -b 100 --bias --fusion -t 2 --rf-stranded >> %s 2>&1' %(genome_dir, kallisto_input_files, results_dir,kallisto_stdout)

    print('Kallisto run command:%s' %kallisto_cmd)
    #subprocess.call(kallisto_cmd, stdout=False)
    os.system(kallisto_cmd)
    
 ####################### Create Kallisto index ###############################
def kallisto_index(org):
    print
    print('\033[33mRunning kallisto index! \033[0m')
    # Organism specific kalisto indexes
    if org == "3H11":
        kallistoindex_cmd = 'kallisto index -i %s/3H11_transcriptome_kallistoindex %s' %(genome_dir,transcriptome_file)
    if org == "R12":
        kallistoindex_cmd = 'kallisto index -i %s/R12_transcriptome_kallistoindex %s' %(genome_dir,transcriptome_file)    
    print('kallisto index command:%s' %(kallistoindex_cmd))
    os.system(kallistoindex_cmd)
   

####################### Running the Pipeline ###############################    
def run_pipeline():
    folder_count = 1
    #data_folders = get_data()

    # Loop through each data folder
    #for data_folder in data_folders:
    folder_name = data_folder.split('/')[-1]
    print
    print
    print('\033[33mProcessing Folder: %s\033[0m' %(folder_name))

    # Get the list of first file names in paired end sequences
    first_pair_files = glob.glob('%s/*_R1*.fastq*' %(data_folder))
    #second_pair_files = glob.glob('%s/_R2*.fastq*' %(data_folder))

    # Program specific results directories
    data_trimmed_dir = "%s/%s/trimmed" %(data_dir,folder_name)
    fastqc_dir = "%s/%s/fastqc_results" %(data_dir,folder_name)
    results_dir = "%s/rnaseq_results_kallisto/%s/%s" %(run_dir,org,folder_name)
    log_file = "%s/%s_run_log.txt" %(results_dir,folder_name)
    kallisto_stdout = "%s/%s_kallisto_stdout.txt" %(results_dir,folder_name)
   

    
    # Run create directories function to create directory structure
    create_dirs(data_trimmed_dir,fastqc_dir,results_dir)
 

    # Loop through each file and create filenames
    file_count = 1
    for first_pair_file in first_pair_files:
        first_file_name_full = first_pair_file.split('/')[-1]
        
        second_pair_file = first_pair_file.replace('_R1', '_R2')
        second_file_name_full = second_pair_file.split('/')[-1]
        file_ext = first_pair_file.split('.')[-1]

        print ('\033[32m Processing File: %s of %s (%s)\033[0m' %(file_count, len(first_pair_files), first_file_name_full ))

        first_file_name = re.split('.fastq|.fastq.gz',first_file_name_full)[0]
        second_file_name = re.split('.fastq|.fastq.gz',second_file_name_full)[0]
        print('first_file_name:%s, second_file_name:%s' %(first_file_name,second_file_name))

        # Collect Sample attributes
        exp_name = folder_name
        print("exp_name: %s" %(exp_name))
        lane = first_file_name.split("_")[-3]
        print("Lane: %s" %(lane))
        sample_id = re.split('.fastq|.fastq.gz', first_file_name)[0]
        print("sample_id: %s" %(sample_id))
        
        
        # 01. Run TrimGalore
        trimgalore(first_pair_file,second_pair_file,folder_name,sample_id,data_trimmed_dir,fastqc_dir)

        file_count = file_count + 1
        #sys.exit()
        
    # Collect trimmed data
    first_pair_group,second_pair_group,kallisto_input_files = collect_trimmed_data(data_trimmed_dir)
    
    # Run salmon if you prefer over to kallisto
    #run_salmon(first_pair_group,second_pair_group,results_dir)
    #sys.stdout = open(kallisto_stdout, 'w') 
    # Run kallisto
    run_kallisto(first_pair_group,second_pair_group,results_dir,kallisto_input_files,org,kallisto_stdout)
    
    
    folder_count = folder_count + 1
    #sys.exit()
    return data_trimmed_dir,fastqc_dir,results_dir,log_file,kallisto_stdout


#salmon_index()
# Create index files if they dont exist.
if not os.path.exists('%s/%s_transcriptome_kallistoindex' %(genome_dir,org)):
    kallisto_index(org) # just need to do it once.
else:
    print ('\033[32m Transcriptome index exists, not creating \033[0m') 


data_trimmed_dir,fastqc_dir,results_dir,log_file,kallisto_stdout = run_pipeline()
   


# with open(log_file,'w') as f:
#     f.write("\n Sample name: ")
#     f.write(folder_name)
#     f.write("\n Pipeline version: %s" %pipeline_version)
#     f.write("\n Transcriptome Reference Files: ")
#     f.write(transcriptome_file)
#     f.write("\n Results directory: ")
#     f.write(results_dir)
#     f.write("\n Data directory: ")
#     f.write(data_dir)
#     f.flush()
#     f.write("\n ------------------------------\n")

#     f.write("\n Kallisto version:\n")
#     f.flush()
#     subprocess.call(['kallisto', 'version'], stdout=f)
#     f.write("\n ------------------------------\n")
#     print()

#     f.write("\n Trim Galore version:\n")
#     f.flush()
#     subprocess.call(['trim_galore', '--version'],stdout=f)
#     f.write("\n ------------------------------\n")
#     print()



