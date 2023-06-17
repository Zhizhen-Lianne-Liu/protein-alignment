import numpy as np

with open('blosum.txt', 'r') as f:
    blosum = [[num for num in line.split()] for line in f]

u = -12

def edit_distance(s, t):
    I = len(s) + 1
    J = len(t) + 1
    D = [[0 for i in range(J)] for j in range(I)]
    path = [[[0, 0] for i in range(J)] for j in range(I)]
    
    for i in range (I):
        D[i][0] = u

    for j in range(J):
        D[0][j] = u

    for i in range(1, I):
        for j in range(1, J):
            m1 = D[i-1][j] -8
            m2 = D[i][j-1] -8
            m3 = D[i-1][j-1] + score(s[i-1], t[j-1])
            D[i][j] = max([m1, m2, m3])
            index = np.argmax([m1, m2, m3])

            if index == 0:
                path[i][j] = [i-1, j]
            elif index == 1:
                path[i][j] = [i, j-1]
            else:
                path[i][j] = [i-1, j-1]

#-----------ignore -------------------------------
            # if index == 0:
            #     edit_seq = edit_seq + 'D'
            #     salign = salign + s[i-1]
            #     talign = talign + ' '

            # elif index == 1:
            #     edit_seq = edit_seq + 'I'
            #     salign = salign + ' '
            #     talign = talign + t[i-1]
            # else:
            #     if compare(s[i-1], t[j-1]) == 0:
            #         edit_seq  = edit_seq + 'M'
            #     else:
            #         edit_seq = edit_seq + 'R'
            #     salign = salign + s[i-1]
            #     talign = talign + t[i-1]
#-------------------------------------------------------
            
    
    align_path(path, s, t)
    return D[len(s)][len(t)], D, path

def align_path(P, s, t):
    edit_seq = ''
    salign = ''
    talign = ''
    current = [len(s), len(t)]
    count = 0
    while current != [0, 0] and count < 50:
        count += 1
        if P[current[0]][current[1]] == [current[0]-1, current[1]-1]:
            if s[current[0]-1] == t[current[1]-1]:
                edit_seq = 'M' + edit_seq
                
            else:
                edit_seq = 'R' + edit_seq
            salign = s[current[0]-1] + salign
            talign = t[current[1]-1] + talign
            current = list(np.subtract(np.array(current), np.array([1, 1])))
        elif P[current[0]][current[1]] == [current[0]-1, current[1]]:
            edit_seq = 'D' + edit_seq
            salign = s[current[0]-1] + salign
            talign = ' ' + talign
            current = list(np.subtract(np.array(current), np.array([1, 0])))
        else:
            edit_seq = 'I' + edit_seq
            salign = ' ' + salign
            talign = t[current[1]-1] + talign
            current = list(np.subtract(np.array(current), np.array([0, 1])))
    print(edit_seq)
    print(salign)
    print(talign)



def score(a, b):
    indexa = blosum[0].index(str(a))
    indexb = blosum[0].index(str(b))
    return int(blosum[indexa][indexb])

def main():
    # some changes were made to protein and blossum files to make reading them easier
    # contrary to before, now the higher the score, the better the alignment.
    file = open('proteins.txt', "r")
    content = file.readlines()
    content = [line.rstrip() for line in content]
    file.close()

    
    
    #print(content)

    proteinA = content[1]
    proteinB = content[3]

    v, D, P = edit_distance(proteinA, proteinB)
    

    print()

if __name__ == '__main__':
    main()

