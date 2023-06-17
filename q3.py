import numpy as np

def edit_distance(s, t):
    I = len(s) + 1
    J = len(t) + 1
    D = [[0 for i in range(J)] for j in range(I)]
    path = [[[0, 0] for i in range(J)] for j in range(I)]
    
    for i in range (I):
        D[i][0] = i

    for j in range(J):
        D[0][j] = j

    for i in range(1, I):
        for j in range(1, J):
            m1 = D[i-1][j] + 1
            m2 = D[i][j-1] + 1
            m3 = D[i-1][j-1] + compare(s[i-1], t[j-1])
            D[i][j] = min([m1, m2, m3])
            index = np.argmin([m1, m2, m3])

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
    while current != [0, 0] and count < 500:
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
    print(edit_seq[:50])
    print(salign[:50])
    print(talign[:50])



def compare(a, b):
    if a == b:
        return 0
    else:
        return 1

def main():
    file = open('proteins.txt', "r")
    content = file.readlines()
    content = [line.rstrip() for line in content]
    #print(content)

    proteinA = content[1]
    proteinB = content[3]

    ed, D, P = edit_distance(proteinA, proteinB)
    #ed, D, P = edit_distance('shesells', 'seashellsa')

    print()

if __name__ == '__main__':
    main()