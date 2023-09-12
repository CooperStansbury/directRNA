import pandas as pd
import numpy as np
import pyranges as pr
import sys


if __name__ == "__main__":
    gtfPath = sys.argv[1]
    outPath = sys.argv[2]

    # load the gtf as a table
    gr = pr.read_gtf(gtfPath)
    gf = gr.as_df()
    
    gf.to_csv(outPath, index=False)
    
    
    
    
 