'''
Module 'check_dependencies' - function to check if all the program dependencies are satisfied 

 A) Python modules
 B) External Programs
 C) files

'''

## DEPENDENCIES ##
# External
import pkg_resources
import shutil
import os


## FUNCTIONS ##


def missing_python_dependencies():
    '''
    check if all the python dependencies are satisfied 
    if not print what is missing and how can be obtained
    
    Input:
        
    Output: true / false
        
    ''' 
    ## 1 list of needed python packages
    needed_python_packages=["numpy", "scipy", "pysam", "mappy", "cigar", "pybedtools"]
    
    ## 2. load packages information
    installed_packages = []   
    dists = [str(d).split(" ")[0] for d in pkg_resources.working_set]

    ## 3. create a list 
    for i in dists:
        installed_packages.append(i)
      
    ## 4. Select not founded packages
    missing_packages=set(needed_python_packages)-set(installed_packages)

    ## 5. check if the set is not empty 
    if bool(missing_packages):
        print("\n\n\n*** ERROR **** Meiga dependencies not satisfied\nPython packages missing : install them with")
        
        ## 6. print list of shell commands to install the missing packages
        for missing_package in missing_packages:
            print("pip3 install "+ missing_package)
        ## 7. exit until fixed
        print("\n\n")
        return True

    ## all 
    return False

def missing_program_dependencies(): 
    '''
    check if all the program dependencies are satisfied 
    if not print what is missing and how can be obtained
    
    Input:
        
    Output: true / false
        
    '''
    missing_external_programs="\n\n\n*** ERROR *** Meiga needs the folowing programs\n"
    is_missing_programs=False

     ## 1 list of needed python packages
    needed_external_programs=["samtools"] #,"bedtools","perl","annovar","bwa-mem","racoon","mimimap2","ncbi-blast"]

    for needed_external_program in needed_external_programs:
        if shutil.which(needed_external_program)== None:
            missing_external_programs+=needed_external_program+"\n"
            is_missing_programs=True
    
    if is_missing_programs:
        print(missing_external_programs)
        return True

    ## all in order
    return False

def missing_needed_files(list_paths):
    '''
    check if all the program dependencies are satisfied 
    if not print what is missing and how can be obtained
    
    Input:
        
    Output: true / false
        
    '''
    missing_files_directories="\n\n\n*** ERROR *** Meiga can't locate the folowing files/directories\n"
    is_missing_files_directories=False

    for path in list_paths:
        if path!=None and not os.path.exists(path):
            missing_files_directories+=path+"\n"
            is_missing_files_directories=True

    #special search for retrotransposons_repeatMasker.bed
    #TODO unify file names 

    if not os.path.exists(list_paths[2]+"/retrotransposons_repeatMasker.bed"):
         missing_files_directories+=list_paths[2]+"/retrotransposons_repeatMasker.bed\n"
         is_missing_files_directories=True

    
    if is_missing_files_directories:
        print(missing_files_directories)
        return True

    ## all in order
    return False

