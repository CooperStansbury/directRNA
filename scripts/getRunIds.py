import pandas as pd
import os


def getRids(dirPath):
    """A function to return the run Ids and sample types of 
    each file in an input directory """
    runIds = []
    
    for f in os.listdir(dirPath):
        fName = f.replace(".fastq.gz", "")
        runIds.append(fName)
        
    return runIds
    