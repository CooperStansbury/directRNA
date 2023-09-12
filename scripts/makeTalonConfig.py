import pandas as pd
import numpy as np
import sys


if __name__ == "__main__":
    platform = sys.argv[1]
    sample = sys.argv[2]
    outpath = sys.argv[3]
    fileList = sys.argv[4:]

    # Dataset config file: dataset name, sample description, platform, sam file (comma-delimited)
    new_rows = []
    for file in fileList:
        fileId = file.split("/")[-1].replace("_labeled.sam", "")
        row = {
            'name' : fileId,
            'sample' : sample,
            'platform' : platform,
            'path' : file,
        }
        new_rows.append(row)

    cf = pd.DataFrame(new_rows)

    # write the config file
    cf.to_csv(outpath, index=False, header=False)
    
    
    
 