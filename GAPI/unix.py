'''
Module 'unix' - Contains wrappers to unix commands
'''

## DEPENDENCIES ##
# External
import os
import subprocess

# Internal
from GAPI import log

## FUNCTIONS ##

def mkdir(path):
    '''
    Create directory

    Note: improve function to be able to create lists of directories

    Input:
        1. path: directory to be created
    '''

    exist = os.path.isdir(path)

    # Only attempt to create directory if it does not exists
    if not exist: 
        try:  
            command = ['mkdir', '-p', path]
            subprocess.call(command)
            
        except OSError:  
            step = 'ERROR'
            msg = "Creation of the directory %s failed" % path
            log.step(step, msg)


def rm(paths):
    '''
    Delete set of files/directories. Directories are deleted recursively

    Input:
        1. files: list containing file/directory paths to be deleted
    '''

    for path in paths:

        exist = os.path.exists(path)

        # Only attempt to delete file if it does exists
        if exist: 
            try:  
                command = ['rm', '-r', path]
                subprocess.call(command)

            except OSError:  
                step = 'ERROR'
                msg = "Deletion of the file %s failed" % filePath
                log.step(step, msg)