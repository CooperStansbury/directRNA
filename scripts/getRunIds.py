import pandas as pd
import os


def getRids(dirPath):
    """A function to return the run Ids and sample types of 
    each file in an input directory """
    sTypes = []
    runIds = []
    
    for f in os.listdir(dirPath):
        sampleType = f.split("_")[0]
        fName = f.replace(".fastq.gz", "").replace("poly_", "").replace("total_", "")
        sTypes.append(sampleType)
        runIds.append(fName)
        
    return sTypes, runIds
    