#!/usr/bin/python

import click
import os
import numpy as np
import pickle
import pandas as pd

@click.command()
@click.option("--pfile","-pf",help="feature file in special format for predicting")
@click.option("--region","-r",help="functional region")
@click.option("--outp","-o",help="outprefix")

def main(pfile, region, outp):
    # read data file
    print("Loading data...")
    Data = pd.read_csv(pfile, encoding='gbk', sep='\t')
    Features = Data[Data.columns[2:len(Data.columns)]]
    Sample = Data['Var']
    size = len(Data)
    print("Data is ready!\nSample size is "+str(size))

    # predicting
    print("Loading model...")
    pwdir = os.path.dirname(os.path.abspath(__file__))
    if '5' in region:
        pklfile = pwdir+'/model/GDBoost_for_5UTR.pickle'
    elif '3' in region:
        pklfile = pwdir+'/model/GDBoost_for_3UTR.pickle'
    model = pickle.load(file=open(pklfile, 'rb'))
    print(model)
    print("Model is ready!")
    print("Start predicting...")
    proba = model.predict_proba(Features)
    
    # output
    outfile = outp+'_pred.txt'
    print("Prediction Result are written into "+outfile)
    out = open(outfile, "w")
    out.write("Sample\tProb\n")
    for i in range(0, len(Sample)):
        oinf = str(Sample[i])+'\t'+str(proba[i][1])
        out.write(oinf+'\n')
    out.close()
    
if __name__ == "__main__":
    main()
