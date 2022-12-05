import pandas as pd
import os


def getRids(dirPath):
    """A function to return the run Ids and sample types of 
    each file in an input directory """
    sTypes = []
    runIds = []
    
    for f in os.listdir(dirPath):
        sampleType = f.split("_")[0]
        fName = "".join(f.split("_")[1:]).replace(".fastq.gz", "")

        sTypes.append(sampleType)
        runIds.append(fName)
        
    return sTypes, runIds
    