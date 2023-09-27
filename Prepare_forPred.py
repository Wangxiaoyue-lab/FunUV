#!/usr/bin/python
# coding=UTF-8

import click
import re
import os
import math
import numpy as np
import gzip

@click.command()
@click.option("--vfile","-v",help="variation file in special format")
@click.option("--region","-r",help="functional region")

def main(vfile, region):
    if '5' in region:
        idx = [0,1,2,3,4,5,6,7,8,9,10,30,31,34,35,38,39,43,50,51,55,56]
    elif '3' in region:
        idx = [0,1,2,3,4,5,6,14,15,25,26,29,30,42,46,47,49,61]
    with open(vfile) as vf:
        for vline in vf:
            vinf = vline.strip().split('\t')
            out = []
            for i in range(0, len(vinf)):
                if i in idx:
                    out.append(vinf[i])
            print('\t'.join(out))

if __name__ == "__main__":
    main()
