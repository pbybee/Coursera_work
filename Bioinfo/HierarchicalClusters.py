import math
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import copy
import numpy as np
import itertools


def EucDist(v,w,m):
    # sum([(vi-wi)**2 for vi in v for wi in w])
    dist = math.sqrt(sum([(v[i]-w[i])**2 for i in range(m)]))
    return dist


if __name__ == "__main__":
    #read in data
    filename = os.path.dirname(__file__) + "/HierData.txt"
    with open(filename, "r") as InputData:
        distance = int(InputData.readline().strip())
        data = {}
        for idx,item in enumerate(InputData.read().splitlines()):
            data[idx] = list(itertools.takewhile(lambda x: x>0.0, map(float,item.split(' '))))

    n = 7
    clusters = list(range(n))
    depth = {}
    while len(clusters)>1:
        minValue = 0
        for values in data.values():
            if values:
                if minValue<max(values):
                    minValue = max(values)
        for key,value in data.items():
            if value:
                if min(data[key]) < minValue:
                    minValue = min(data[key])
                    minIdx = data[key].index(minValue)
                    minKey = key

        #merge clusters
        depth['-'.join([str(minKey),str(minIdx)])] = minValue/2
        data['-'.join([str(minKey),str(minIdx)])] = 0
        del data[minKey]
        del data[minIdx]
        for key,value in data.items():
            i = 1


    for item in nodes[0].points:
        print(item+'\n')

