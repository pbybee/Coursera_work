import math
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import copy
import numpy as np

#Soft clustering assigns data points to centers with a weight or responsibility.
def EucDist(v,w,m):
    # sum([(vi-wi)**2 for vi in v for wi in w])
    dist = math.sqrt(sum([(v[i]-w[i])**2 for i in range(m)]))
    return dist


if __name__ == "__main__":
    #read in data
    filename = os.path.dirname(__file__) + "/Soft_KMeansData.txt"
    with open(filename, "r") as InputData:
        k,m = map(int,InputData.readline().strip().split(' '))
        beta = float(InputData.readline().strip())
        data = []
        for item in InputData.read().splitlines():
            data.append(list(map(float,item.split(' '))))

    centers = []
    currentCenters = np.array(data[0:k])
    clusters = defaultdict(list)
    hiddenMatrix = np.zeros((k,len(data)))

    while not np.array_equal(np.round(centers,3),np.round(currentCenters,3)):

        centers = copy.deepcopy(currentCenters)
        clusters = defaultdict(list)

        for j in range(len(data)):
            secondTerm = 0
            for center in centers:
                secondTerm += math.exp(-beta*EucDist(data[j],center,m))
            for i in range(len(centers)):
                hiddenMatrix[i][j] = math.exp(-beta*EucDist(data[j],centers[i],m))/secondTerm


        for i in range(len(centers)):
            for j in range(m):
                currentCenters[i][j] = sum([HM*dat[j] for HM,dat in zip(hiddenMatrix[i],data)])/sum(hiddenMatrix[i])



    for item in centers:
        print(' '.join(['%.3f' % elem for elem in item]))
