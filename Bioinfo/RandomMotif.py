__author__ = 'ptb-t440'

import pandas as pd
import os
from random import randint
from MotifAlgorithms import MotifAlgorithms

#Not used in this exercise. Lead up to Gibbs sampling
def RandomMotifSearch(nucStrings, k, t):
    randomInts = [randint(0,len(nucStrings[0])-k) for a in range(t)]
    motifs = [nucStrings[i][r:r+k] for i,r in enumerate(randomInts)]

    bestMotifs = [MotifAlgorithms.KmerScore(motifs), motifs]
    while 1:
        currentProfile = MotifAlgorithms.ProfilPseudocounts(motifs)
        motifs = list(MotifAlgorithms.MostProbableKmers(dna,k,currentProfile) for dna in nucStrings)
        currentScore = MotifAlgorithms.KmerScore(motifs)
        if currentScore < bestMotifs[0]:
            bestMotifs = [currentScore,motifs]
        else:
            return (bestMotifs)


def GibbsSampler(dna,k,t,N):
    # Randomly generate k-mers from each sequence in the dna list.
    rand_ints = [randint(0,len(dna[0])-k) for a in range(t)]
    motifs = [dna[i][r:r+k] for i,r in enumerate(rand_ints)]

    # Initialize the best score as a score higher than the highest possible score.
    best_score = [MotifAlgorithms.KmerScore(motifs), motifs]

    # Iterate motifs.
    for i in range(N):
        r = randint(0,t-1)
        current_profile = MotifAlgorithms.ProfilPseudocounts([motif for index, motif in enumerate(motifs) if index!=r])
        # print 'a: ', motifs
        motifs = [MotifAlgorithms.MostProbableKmers(dna[index],k,current_profile) if index == r else motif for index,motif in enumerate(motifs)]
        # print 'b: ', motifs
        current_score = MotifAlgorithms.KmerScore(motifs)
        if current_score < best_score[0]:
            best_score = [current_score, motifs]

    return best_score


def main():

    #Read in the data and apply the Gibbs sampler.
    filename = os.path.dirname(__file__) + "/RandData.txt"
    with open(filename, "r") as randData:
        #k is kmer nucleotides, t is number of strings, N is number of cycles
        k,t,N = map(int, randData.readline().split(' '))
        nucStrings = list(randData.read().splitlines())
        bestMotifs = [k*t, None]
        for i in range(20):
            newbestMotifs = GibbsSampler(nucStrings, k, t, N)

            if newbestMotifs[0] < bestMotifs[0]:
                bestMotifs = newbestMotifs
        print (bestMotifs[0])
        for item in bestMotifs[1]:
            print (item+'\n')





if __name__ == "__main__":
    main()

