import pandas as pd
import os
import numpy as np
from collections import defaultdict, OrderedDict


#This file find the Eurler path through a graph
def main ():
    outFile = open("euler_Output.txt", 'w+')

    filename = os.path.dirname(__file__) + "/eulerData.txt"

    with open(filename, "r") as nucData:
        nucString = nucData.read().splitlines()

        EulerGraph = defaultdict(list)
        LenDict = defaultdict(int)
        for word in range(len(nucString)):
            tempItem = nucString[word].split(' -> ')
            EulerGraph[int(tempItem[0])] = [int(x) for x in tempItem[1].split(',')]
        for key,value in EulerGraph.items():
            LenDict[key] = len(value)

        EulerPath = list()

        #If there are multiple indexes for the same item value, find all indexes and return in t a list
        multipleIdx = lambda searchList, item: [idx for idx in range(len(searchList)) if searchList[idx]==item]

        #Do it the first time
        key, value = EulerGraph.popitem()
        for v in value:
            EulerPath.extend([key,v])

        #while there are still keys in the dictionary
        while (len(EulerGraph.keys())):
            key, value = EulerGraph.popitem()
            valTimes = [j for j in value if j in EulerPath]
            if key in EulerPath:
                keyIndices = multipleIdx(EulerPath, key)
                for idx in keyIndices:
                    if not LenDict.get(value[-1]):
                        #then it is the last element
                        EulerPath.append(value.pop())
                    elif EulerPath.count(value[-1]) < LenDict[value[-1]]:
                        EulerPath.insert(idx+1, value.pop())
                    else:
                        value.pop()

            elif valTimes:
                for j in valTimes:
                    valIndices = multipleIdx(EulerPath,j)
                    for vidx in valIndices:
                        if EulerPath.count(j) < LenDict[j] or key not in EulerPath:
                            EulerPath.insert(vidx,key)

            while value:
                if EulerPath.count(value[-1]) < LenDict[value[-1]]:
                    EulerPath.extend([key, value.pop()])
                else:
                    value.pop()

        # outFile.write('->'.join(map(str,EulerPath)) + '\n')
        print('->'.join(map(str,EulerPath)) + '\n')




if __name__ == "__main__":
    main()
