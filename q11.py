
import numpy as np

with open('blosum.txt', 'r') as f:
    blosum = [[num for num in line.split()] for line in f]

def score(a, b):
    indexa = blosum[0].index(str(a))
    indexb = blosum[0].index(str(b))
    return int(blosum[indexa][indexb])


def V_sub(s, t):
    V_sfx = [[0 for j in range(len(t)+1)]for i in range(len(s)+1)]

    for i in range(1, len(s) + 1):
        for j in range(1, len(t) + 1):
            V_sfx[i][j] = max(0, V_sfx[i-1][j-1] + score(s[i-1], t[j-1]), V_sfx[i-1][j] -2, V_sfx[i][j-1]-2)

    return max(max(V_sfx))



def main():
    # some changes were made to protein and blossum files to make reading them easier
    # contrary to before, now the higher the score, the better the alignment.
    file = open('proteins.txt', "r")
    content = file.readlines()
    content = [line.rstrip() for line in content]
    file.close()
    
    #print(content)

    proteinC = content[5]
    proteinD = content[7]
    
    Vsub = V_sub(proteinC, proteinD)
    print(Vsub)

    print()

if __name__ == '__main__':
    main()