import pandas as pd
import os
import numpy as np

if __name__ == "__main__":
    filename = os.path.dirname(__file__)+"/AlignData.txt"
    bloFile = os.path.dirname(__file__)+"/PAM250_Dict.txt"
    with open(bloFile) as bloData:
        blosum = eval(bloData.read())

    sigma = 5
    with open(filename) as PathData:
        v,w = [line.strip() for line in PathData.readlines()]

    S = np.zeros((len(v)+1,len(w)+1), dtype=int)
    backTrack = np.zeros((len(v)+1,len(w)+1), dtype=int)

    for i in range(1, len(v)+1):
        S[i][0] = -i*sigma
    for j in range(1, len(w)+1):
        S[0][j] = -j*sigma

    #apply blosum penalty
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            scores = [S[i-1][j] - sigma, S[i][j-1] - sigma, S[i-1][j-1] + blosum[v[i-1]][w[j-1]]]
            S[i][j] = max(scores)
            backTrack[i][j] = scores.index(S[i][j])

    # Quick lambda function to insert indels.
    insert_indel = lambda word, i: word[:i] + '-' + word[i:]
    # Initialize the aligned strings as the input strings.
    v_aligned, w_aligned = v, w

    # Get the position of the highest scoring cell in the matrix and the high score.
    i, j = len(v), len(w)
    max_score = str(S[i][j])

    # Get the position of the highest scoring cell in the matrix and the high score.
    # Initialize the aligned strings as the input strings up to the position of the high score.
    while i*j != 0:
        if backTrack[i][j] == 0:
            i -= 1
            w_aligned = insert_indel(w_aligned, j)
        elif backTrack[i][j] == 1:
            j -= 1
            v_aligned = insert_indel(v_aligned, i)
        else:
            i -= 1
            j -= 1

    # Prepend the necessary preceeding indels to get to (0,0).
    for repeat in range(i):
        w_aligned = insert_indel(w_aligned, 0)
    for repeat in range(j):
        v_aligned = insert_indel(v_aligned, 0)

    editDistance = 0
    for i in range(len(w_aligned)):
        if w_aligned[i] != v_aligned[i]:
            editDistance += 1

    print(str(editDistance) + '\n' + v_aligned + '\n' + w_aligned)