import pandas as pd
import os
import numpy as np
from collections import defaultdict
from MotifAlgorithms import MotifAlgorithms

#map the kmers to a graph with a list of the kmers and a set or dictionary of the kmers prefix suffix overlaps

#make dictionary with key as the kmer, and value as the intersections or overlaps

def main ():
    # outFile = open("/Mismatch_Output.txt", 'w+')

    filename = os.path.dirname(__file__) + "/GenomePathData.txt"
#     filename = "/GenomePathData.txt"

    with open(filename, "r") as nucData:
        nucString = nucData.read().splitlines()

        prefixgraph = MotifAlgorithms.KmerGraph(nucString)
        for item in prefixgraph:
            print(item[0] + " -> " + item[1])


if __name__ == '__main__':
    main()

