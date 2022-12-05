import pandas as pd
import numpy as np
import sys

def getGenes(gtfPath):
    """A function to return the gene names by parsing the 
    metadata in the GTF file
    
    args:
        : gtfPath (str): full path to the GTF file

    returns:
        : df (pd.DataFrame): the ensembleID and the gene name 
        for every entry in the GTF file
    """
    
    with open(gtfPath) as f:
        gtf = list(f)
    
    # exclude the comment lines
    gtf = [x for x in gtf if not x.startswith("#")]
    
    newRows = []
    
    
    for line in gtf:
        if 'gene_id' and 'gene_name' in line:
            meta = line.split("\t")[-1]
            geneId = meta.split('gene_id')[1].split(";")[0].replace('"', "").strip()
            geneName = meta.split('gene_name')[1].split(";")[0].replace('"', "").strip()
            
            geneBiotype = meta.split('gene_biotype')[1].split(";")[0].replace('"', "").strip()
            
            if 'transcript_id' in line:
                transcriptId = meta.split('transcript_id')[1].split(";")[0].replace('"', "").strip()
                transcriptName = meta.split('transcript_name')[1].split(";")[0].replace('"', "").strip()
            else:
                transcriptId = 'None'
                transcriptName = 'None'

            row = {
                'geneId' : geneId,
                'geneName' : geneName,
                'transcriptId' : transcriptId,
                'transcriptName' : transcriptName,
                'geneBiotype' : geneBiotype
            }
            newRows.append(row)
            
    df = pd.DataFrame(newRows)
    return df
    


if __name__ == "__main__":
    gtfPath = sys.argv[1]
    outPath = sys.argv[2]
    
    df = getGenes(gtfPath)
    df = df.drop_duplicates()
    df.to_csv(outPath, index=False)
    
    
    
    
 