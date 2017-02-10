import pandas as pd
import os
import numpy as np
from collections import defaultdict, OrderedDict
import operator
from itertools import combinations, product
from functools import reduce
from math import sqrt
from random import randint


class MotifAlgorithms:

    # def __init__ (self):
    #     pass

    @staticmethod
    def MedianString(nucString, k):
        distDict = defaultdict(int)
        for kmer in product("ACGT",repeat=k):
            kmer = ''.join(kmer)
            kmerScore = list()
            for kmerList in list(MotifAlgorithms.Composition(DNA) for DNA in nucString):
                kmerScore.append(min([MotifAlgorithms.HammingDistance(kmer, motif) for motif in kmerList]))
            distDict[kmer] = sum(kmerScore)

        dictTuple = sorted(distDict.items(), key=lambda x: (x[1],x[0]), reverse=True)

        return dictTuple

    @staticmethod
    def GreedyMotif(nucStrings, k, t):
        best_score = [t*k, None]
        for i in range(len(nucStrings[0])-k+1):
            # Initialize the motifs as each k-mer from the first dna sequence.
            motifs = [nucStrings[0][i:i+k]]
            current_profile = MotifAlgorithms.Profile(motifs)

            # Find the most probable k-mer in the next string.
            for j in range(1,t):
                motifs.append(MotifAlgorithms.MostProbableKmers(nucStrings[j],k,current_profile))
                current_profile = MotifAlgorithms.Profile(motifs)

            # Check to see if we have a new best scoring list of motifs.
            current_score = MotifAlgorithms.KmerScore(motifs)
            if current_score < best_score[0]:
                best_score = [current_score, motifs]
        return (best_score[1])

    @staticmethod
    def MostProbableMotif(nucString, k, profileProb):
        kmers = list(kmer for kmer in [nucString[idx:idx+k] for idx in range(len(nucString)-k+1)])
        kmerDict = defaultdict(float)
        maxProb = 0
        for kmer in kmers:
            kmerProb = list()
            for idx in range(len(kmer)):
                kmerProb.append(profileProb[kmer[idx]][idx])
            kmerDict[kmer] = reduce(operator.mul, kmerProb)
            if kmerDict[kmer] > maxProb:
                maxProb = kmerDict[kmer]
                maxKmer = kmer
        return (maxKmer)

#############################Before Week 3############################

    @staticmethod
    def ReverseComposition (nucString):

        mergedString = ""
        for idx, item in enumerate(nucString):
            if idx == 0:
                mergedString = item
            else:
                mergedString += item[-1]

        return (mergedString)

    @staticmethod
    def DeBruijnGraph (nucString, kmer):
        compList = MotifAlgorithms.Composition(nucString,kmer)
        # compList = nucString
        sortList = []
        debruGraph= defaultdict(list)
        for idx in range(len(compList)):
            sortList.append(MotifAlgorithms.PatternToNumber(compList[idx]))
        sortList.sort()
        sortedKmers = []
        for idx in range(len(sortList)):
            sortedKmers.append(MotifAlgorithms.NumberToPattern(sortList[idx], kmer))

            #find the possible suffixes for the current kmer
            Overlap = MotifAlgorithms.PatternToNumber(sortedKmers[-1][1:])
            matchIdx = 0
            while matchIdx <= 3:
                matchNum = 4*Overlap+matchIdx
                if (matchNum in sortList  and sortList.index(matchNum) != idx):
                    matchStr = MotifAlgorithms.NumberToPattern(matchNum, kmer)
                    debruGraph[sortedKmers[idx][:kmer-1]].append(MotifAlgorithms.NumberToPattern(matchNum, kmer)[:kmer-1])
                matchIdx +=1

        return debruGraph

    @staticmethod
    def KmerGraph (nucString):
        sortList = []
        adjacentList= []
        kmer = len(nucString[0])
        for idx in range(len(nucString)):
            sortList.append(MotifAlgorithms.PatternToNumber(nucString[idx]))
        sortList.sort()
        sortedKmers = []
        for idx in range(len(sortList)):
            sortedKmers.append(MotifAlgorithms.NumberToPattern(sortList[idx], kmer))

            #find the possible suffixes for the current kmer
            Overlap = MotifAlgorithms.PatternToNumber(sortedKmers[-1][1:])
            matchIdx = 0
            while matchIdx <= 3:
                if (4*Overlap+matchIdx in sortList):
                    adjacentList.append([sortedKmers[idx],MotifAlgorithms.NumberToPattern(4*Overlap+matchIdx, kmer)])
                matchIdx +=1
        return adjacentList

    @staticmethod
    def Composition (nucString, kmer):

        k_mers = []
        for idx in range(0, len(nucString)-kmer+1):
            k_mers.append(nucString[idx:idx+kmer])


        return k_mers

    @staticmethod
    def MostProbableKmers (dna,k, profileProb):
        # A dictionary relating nucleotides to their position within the profile.
        nuc_loc = {nucleotide:index for index,nucleotide in enumerate('ACGT')}
        # Initialize the maximum probabily.
        max_prob = [-1, None]
        # Compute the probability of the each k-mer, store it if it's currently a maximum.
        for i in range(len(dna)-k+1):
            current_prob = 1
            for j, nucleotide in enumerate(dna[i:i+k]):
                current_prob *= profileProb[j][nuc_loc[nucleotide]]
            if current_prob > max_prob[0]:
                max_prob = [current_prob, dna[i:i+k]]

        return (max_prob[1])

    @staticmethod
    def Profile (kmers):
        prof = []
        for i in range(len(kmers[0])):
            col = ''.join([kmers[j][i] for j in range(len(kmers))])
            prof.append([float(col.count(nuc))/float(len(col)) for nuc in 'ACGT'])
        return (prof)

    @staticmethod
    def ProfilPseudocounts(motifs):
        '''Returns the profile of the dna list motifs.'''
        prof = []
        for i in range(len(motifs[0])):
            col = ''.join([motifs[j][i] for j in range(len(motifs))])
            prof.append([float(col.count(nuc)+1)/float(len(col)+4) for nuc in 'ACGT'])
        return prof



    @staticmethod
    def KmerScore (motifs):
        score = 0
        for i in range(len(motifs[0])):
            motif = ''.join([motifs[j][i] for j in range(len(motifs))])
            score += min([MotifAlgorithms.HammingDistance(motif, homogeneous*len(motif)) for homogeneous in 'ACGT'])
        return score

    #####################SCORING########################
    ###################################################
    @staticmethod
    def HammingDistance(str1, str2):

        if len(str1) != len(str2):
            return 0
        hammingDist = 0
        for idx in range(0,len(str1)):
            if str1[idx] != str2[idx]:
                hammingDist += 1
        return hammingDist

    @staticmethod
    def MismatchAppearances (nucString, searchPattern, d_ham):
        patternLocation = []
        patternLength = len(searchPattern)
        nucStringLength = len(nucString)
        idx = 0
        while (idx < nucStringLength-patternLength+1):
            currentNucString = nucString[idx:(idx+patternLength)]
            hamDist = MotifAlgorithms.HammingDistance(currentNucString, searchPattern)
            if (hamDist <= d_ham):
                patternLocation.append(idx)
            idx += 1

        return patternLocation

    @staticmethod
    def Neighbors(searchPattern, d_ham):
        if d_ham == 0:
            return [searchPattern]
        if len(searchPattern) == 1:
            return ['A','C','G','T']
        Neighborhood = set()
        suffixNeighbor = MotifAlgorithms.Neighbors(searchPattern[1:], d_ham)
        for strItem in suffixNeighbor:
            if MotifAlgorithms.HammingDistance(searchPattern[1:], strItem) < d_ham:
                for nucl in "ACGT":
                    Neighborhood.add(nucl + strItem)
            else:
                Neighborhood.add(searchPattern[0] + strItem)
        return list(Neighborhood)

    @staticmethod
    def FrequentWords_MisMatch (nucString, kmer, d_ham):
        kmer_freq = defaultdict(int)

        for idx in range(0,len(nucString)-kmer+1):
            kmer_freq[nucString[idx:idx+kmer]] += 1
            # neighborhood.append(MotifAlgorithms.Neighbors(nucString[idx:idx+kmer], d_ham))

        mismatchCount = defaultdict(int)
        for k_word, freq in kmer_freq.items():
            neighborhood = MotifAlgorithms.kmer_mismatches(k_word, d_ham)
            for neighbor in neighborhood:
                mismatchCount[neighbor] += freq

        maxCount = max(mismatchCount.values())
        return sorted(kmer for kmer, count in mismatchCount.items() if count == maxCount)

    @staticmethod
    def kmer_mismatches(kmer, d):
        """Returns all k-mers that are within d mismatches of the given k-mer."""
        mismatches = [kmer]  # Initialize mismatches with the k-mer itself (i.e. d=0).
        alt_bases = {'A':'CGT', 'C':'AGT', 'G':'ACT', 'T':'ACG'}
        for dist in range(1, d+1):
            for change_indices in combinations(range(len(kmer)), dist):
                for substitutions in product(*[alt_bases[kmer[i]] for i in change_indices]):
                    new_mistmatch = list(kmer)
                    for idx, sub in zip(change_indices, substitutions):
                        new_mistmatch[idx] = sub
                    mismatches.append(''.join(new_mistmatch))
        return mismatches

    @staticmethod
    def FrequentWords_Anti_Sense (nucString, kmer, d_ham):
        kmer_freq = defaultdict(int)

        for idx in range(0,len(nucString)-kmer+1):
            kmer_freq[nucString[idx:idx+kmer]] += 1
            # neighborhood.append(MotifAlgorithms.Neighbors(nucString[idx:idx+kmer], d_ham))

        mismatchCount = defaultdict(int)
        for k_word, freq in kmer_freq.items():
            neighborhood = MotifAlgorithms.kmer_mismatches(k_word, d_ham)
            for neighbor in neighborhood:
                mismatchCount[neighbor] += freq

        maxCount = 0
        sense_list = []
        antisense_list = []
        for key, value in mismatchCount.items():
            antiKey = MotifAlgorithms.NucReverseComp(key)
            sense_Count = MotifAlgorithms.MismatchAppearances(nucString, key, d_ham)
            antisense_Count = MotifAlgorithms.MismatchAppearances(nucString, antiKey, d_ham)
            curMax = len(sense_Count)+len(antisense_Count)
            if curMax == maxCount:
                sense_list.append(key)
                antisense_list.append(antiKey)
            if curMax > maxCount:
                maxCount = curMax
                sense_list = []
                antisense_list = []
                sense_list.append(key)
                antisense_list.append(antiKey)


        return sense_list, antisense_list


##########################CHAPTER 1 ALGORITHMS################################

    @staticmethod
    def Transcribe(nucString):
        transcribed = ""
        for nucl in nucString:
            if nucl == 'A':
                transcribed += 'A'
            elif nucl == 'T':
                transcribed += 'U'
            elif nucl == 'G':
                transcribed += 'G'
            elif nucl == 'C':
                transcribed += 'C'
        return transcribed

    @staticmethod
    def PatternCount(nucString, searchPattern):
        #self.nucString = nucString

        #define pattern
        #self.searchPattern = searchPattern
        patternLength = len(searchPattern)

        nucStringLength = len(nucString)
        patternCount = 0
        idx = 0
        for idx in range(len(nucStringLength)-patternLength+1):
            currentNucString = nucString[idx:(idx+patternLength)]
            if (currentNucString == searchPattern):
                patternCount += 1
        return patternCount

    @staticmethod
    def ReverseComp (nucString):
        reverseComp = ""
        for nucl in nucString:
            if nucl == 'A':
                reverseComp += 'T'
            elif nucl == 'T':
                reverseComp += 'A'
            elif nucl == 'G':
                reverseComp += 'C'
            elif nucl == 'C':
                reverseComp += 'G'
        """
        ALTERNATIVE METHOD
        """
#            if nucl == 'A':
#                reverseComp = 'T' + reverseComp
#            elif nucl == 'T':
#                reverseComp = 'A' + reverseComp
#            elif nucl == 'G':
#                reverseComp = 'C' + reverseComp
#            elif nucl == 'C':
#                reverseComp = 'G' + reverseComp

        reverseComp = reverseComp[::-1]
        return reverseComp

    @staticmethod
    def PatternIndex (nucString, searchPattern):
        patternLocation = []
        patternLength = len(searchPattern)
        nucStringLength = len(nucString)
        idx = 0
        for idx in range(len(nucStringLength)-patternLength+1):
            currentNucString = nucString[idx:(idx+patternLength)]
            if (currentNucString == searchPattern):
                patternLocation.append(idx)
        return patternLocation

    @staticmethod
    def FrequentWords_naive (nucString, kmer):
        nucStringLength = len(nucString)
        patternCount = 0
        idx = 0
        kmerCount = []

        while (idx <= nucStringLength-kmer):
            currentNucString = nucString[idx:(idx+kmer)]
            kmerCount.append(MotifAlgorithms.PatternCount(nucString, currentNucString))
            idx += 1

        maxKmerCount = max(kmerCount)

        idx = 0
        frequentPatterns = []
        while (idx <= nucStringLength-kmer):
            if (kmerCount[idx] == maxKmerCount):
                if not(nucString[idx:(idx+kmer)] in frequentPatterns):
                    frequentPatterns.append(nucString[idx:(idx+kmer)])
            idx += 1
        return (frequentPatterns)

    @staticmethod
    def ClumpFinding (nucString, kmer, t_thr, L_win):

        FrequentPatterns = []
        Clumps = [0]*(4**kmer)
        subStr = nucString[0:L_win]
        freqArray = MotifAlgorithms.ComputingFreqs(subStr, kmer)

        for idx in range(0,(4**kmer)-1):
            if freqArray[idx] >= t_thr:
                Clumps[idx] = 1

        for idx in range(1, len(nucString) - L_win):
            firstPattern = nucString[idx-1:idx+kmer-1]
            patternIndex = MotifAlgorithms.PatternToNumber(firstPattern)
            freqArray[patternIndex] = freqArray[patternIndex] - 1
            lastPattern = nucString[idx+L_win-kmer:idx+L_win]
            patternIndex = MotifAlgorithms.PatternToNumber(lastPattern)
            freqArray[patternIndex] = freqArray[patternIndex] + 1
            if freqArray[patternIndex] >= t_thr:
                Clumps[patternIndex] = 1

        for clump_idx in range(0,(4**kmer)-1):
            if Clumps[clump_idx] == 1:
                clumpPattern = MotifAlgorithms.NumberToPattern(clump_idx, kmer)
                FrequentPatterns.append(clumpPattern)

        return FrequentPatterns

    @staticmethod
    def FreqArray (nucString, kmer):
        nucStringLength = len(nucString)
        freqArray = [0]*(4**kmer)
        for idx in range(0, nucStringLength-kmer):
            kmerPattern = nucString[idx:idx+kmer]
            patternNumber = MotifAlgorithms.PatternToNumber(kmerPattern)
            freqArray[patternNumber] += 1
        return freqArray

    #Converts a base pair sequence to an alphebetically/ Lexicographically
            #ordered list
    @staticmethod
    def PatternToNumber(nucString):
        if not nucString:
            return 0
        symbol = nucString[-1]
        prefix = nucString [0:-1]

        return (4*MotifAlgorithms.PatternToNumber(prefix) + MotifAlgorithms.SymbolToNumber(symbol))

    @staticmethod
    def NumberToPattern (str_idx, kmer):
        if kmer == 1:
            return MotifAlgorithms.NumberToSymbol(str_idx)
        symbolIndex = str_idx % 4
        symbol = MotifAlgorithms.NumberToSymbol(symbolIndex)
        prefixIndex = (str_idx - symbolIndex) / 4

        prefixPattern = MotifAlgorithms.NumberToPattern(prefixIndex, kmer-1)

        return (prefixPattern + symbol)

    @staticmethod
    def SymbolToNumber(nucl):
        if nucl == 'A':
            return 0
        elif nucl == 'T':
            return 3
        elif nucl == 'G':
            return 2
        elif nucl == 'C':
            return 1

    @staticmethod
    def NumberToSymbol(nucNum):
        if nucNum == 0:
            return'A'
        elif nucNum == 3:
            return'T'
        elif nucNum == 2:
            return'G'
        elif nucNum == 1:
            return'C'

