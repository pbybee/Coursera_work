import math
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import copy


def EucDist(v,w,m):
    # sum([(vi-wi)**2 for vi in v for wi in w])
    dist = math.sqrt(sum([(v[i]-w[i])**2 for i in range(m)]))
    return dist

#Lloyds algorithm finds the center to clusters distance and conversely the clusters to centers distance to determin
#center of gravity for the clusters.
if __name__ == "__main__":
    #read in data
    filename = os.path.dirname(__file__) + "/ClustData.txt"
    with open(filename, "r") as InputData:
        k,m = map(int,InputData.readline().strip().split(' '))
        data = []
        for item in InputData.read().splitlines():
            data.append(list(map(float,item.split(' '))))
    # plt.scatter(*zip(*data))
    # plt.show()

    centers = []
    currentCenters = data[0:k]
    clusters = defaultdict(list)

    while centers != currentCenters:

        centers = copy.deepcopy(currentCenters)
        clusters = defaultdict(list)

        for point in data:
            pointmax = []
            for i in range(len(centers)):
                pointmax.append(EucDist(point,centers[i],m))
            index = pointmax.index(min(pointmax))
            clusters[index].append(point)

        for i, cluster in enumerate(clusters.values()):
            clusterDimensions = zip(*cluster)
            for j in range(m):
                currentCenters[i][j] = sum(clusterDimensions[j])/len(clusterDimensions[j])
        # for i in range(k):
        #     roundCenters = [round(elem,3) for elem in centers[i]]
        #     roundCurrent = [round(elem,3) for elem in currentCenters[i]]
    for item in centers:
        print(' '.join(['%.3f' % elem for elem in item]))

    plt.scatter(*zip(*data))
    plt.scatter(*zip(*centers),color='red')
    plt.show()
