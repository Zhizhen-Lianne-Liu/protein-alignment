#############################
#copy of q6.py with changeable proteins;returns gap alignment score;uses new scoring matrix
#####################################



from cmath import inf
import numpy as np


def score(x, y):
    if x != ' ':
        indexa = blosum[0].index(str(x))
        indexb = blosum[0].index(str(y))
        return int(blosum[indexa][indexb])
    else:
        return 0

# #----one possible alignment------------------
# def create_path(a, b, c, d, e, f, g, s, t):
#     current = [len(s), len(t)]
#     count = 0

#     #stores whether previous step in path is up or left
#     path = [[[0, 0] for i in range(J+1)] for j in range(I+1)]



#     while current != [0, 0] and count < 500:
#         count += 1

#         #for debugging:
#         A = a[current[0]][current[1]]
#         B = b[current[0]][current[1]]
#         C = c[current[0]][current[1]]
#         D = d[current[0]][current[1]]
#         E = e[current[0]][current[1]]
#         F = f[current[0]][current[1]]
#         G = g[current[0]][current[1]]


#         if (a[current[0]][current[1]] == 1 )and ((e[current[0]][current[1]]== 1 and path[current[0]+1][current[1]] == current)):
#             path[current[0]][current[1]] = [current[0]-1, current[1]]
#             current = [current[0]-1, current[1]]
#             while a[current[0]][current[1]] == 1 and e[current[0]][current[1]]== 1:
#                 path[current[0]][current[1]] = [current[0]-1, current[1]]
#                 current = [current[0]-1, current[1]]
#             if d[current[0]+1][current[1]] == 1:
#                 path[current[0]][current[1]] = [current[0]-1, current[1]]
#                 current = [current[0]-1, current[1]]
#             # if d[current[0]][current[1]] == 1:
#             #     path[current[0]-1][current[1]] = [current[0], current[1]]
#             #     path[current[0]-1][current[1]-1] = [0, 0]
#         elif (b[current[0]][current[1]] == 1) and ((g[current[0]][current[1]]== 1 and path[current[0]][current[1]+1] == current)):
#             path[current[0]][current[1]] = [current[0], current[1]-1]
#             current = [current[0], current[1]-1]
#             while b[current[0]][current[1]] == 1 and g[current[0]][current[1]]== 1:
#                 path[current[0]][current[1]] = [current[0], current[1]-1]
#                 current = [current[0], current[1]-1]
#             if f[current[0]][current[1]+1] == 1:
#                 path[current[0]][current[1]] = [current[0], current[1]-1]
#                 current = [current[0], current[1]-1]


#         elif (b[current[0]][current[1]] == 1) and ((g[current[0]][current[1]]== 0 and path[current[0]][current[1]+1] == [0,0])):
#             path[current[0]][current[1]] = [current[0], current[1]-1]
#             current = [current[0], current[1]-1]
#             if f[current[0]][current[1]+1] == 1:
#                 path[current[0]][current[1]] = [current[0], current[1]-1]
#                 current = [current[0], current[1]-1]

#         elif (a[current[0]][current[1]] == 1 ) and ((e[current[0]][current[1]]== 0 and path[current[0]+1][current[1]] == [0,0])):
#             path[current[0]][current[1]] = [current[0]-1, current[1]]
#             current = [current[0]-1, current[1]]
#             if d[current[0]+1][current[1]] == 1:
#                 path[current[0]][current[1]] = [current[0]-1, current[1]]
#                 current = [current[0]-1, current[1]]
#         elif c[current[0]][current[1]] == 1:
#             path[current[0]][current[1]] = [current[0]-1, current[1]-1]
#             current = [current[0]-1, current[1]-1]
        

#     align_path(path, s, t)

  
# def align_path(path,s, t):
#     edit_seq = ''
#     salign = ''
#     talign = ''
#     current = [len(s), len(t)]
#     count = 0
#     while current != [0, 0] and count < 500:

#         count += 1
#         if path[current[0]][current[1]] == [current[0]-1, current[1]-1]:
#             if s[current[0]-1] == t[current[1]-1]:
#                 edit_seq = 'M' + edit_seq
                
#             else:
#                 edit_seq = 'R' + edit_seq
#             salign = s[current[0]-1] + salign
#             talign = t[current[1]-1] + talign
#             current = list(np.subtract(np.array(current), np.array([1, 1])))
#         elif path[current[0]][current[1]] == [current[0]-1, current[1]]:
#             edit_seq = 'D' + edit_seq
#             salign = s[current[0]-1] + salign
#             talign = '-' + talign
#             current = list(np.subtract(np.array(current), np.array([1, 0])))
#         else:
#             edit_seq = 'I' + edit_seq
#             salign = '-' + salign
#             talign = t[current[1]-1] + talign
#             current = list(np.subtract(np.array(current), np.array([0, 1])))

def Find_score(s, t, u):

    global J, I, P, Q, R, a, b, c, d, e, f, g, blosum

    J = len(t)+1
    I = len(s)+1
    P = [[0 for j in range(J)]for i in range(I)]
    R = [[0 for j in range(J)]for i in range(I)]
    Q = [[0 for j in range(J)]for i in range(I)]

    a = [[0 for j in range(J+1)]for i in range(I+1)]
    b = [[0 for j in range(J+1)]for i in range(I+1)]
    c = [[0 for j in range(J+1)]for i in range(I+1)]
    d = [[0 for j in range(J+1)]for i in range(I+1)]
    e = [[0 for j in range(J+1)]for i in range(I+1)]
    f = [[0 for j in range(J+1)]for i in range(I+1)]
    g = [[0 for j in range(J+1)]for i in range(I+1)] #setting all bit arrays to 0

    blosum = [['#', 'a', 'b'], ['a', '1', '-1'], ['b', '-1', '1']]

    



    for j in range(J):
        P[0][j] = -inf
        R[0][j] = u
    for i in range(I):
        Q[i][0] = -inf
        R[i][0] = u

    R[0][0] = 0
    c[I][J] = 1



    for i in range(1, I):
        for j in range(1, J):
            #find min cost of path ending at N[i][j] using edge V[i][j]
            P[i][j] = max([P[i-1][j], R[i-1][j] + u])
            if P[i][j] == P[i-1][j]:
                d[i-1][j] = 1

            if P[i][j] == R[i-1][j] + u:
                e[i-1][j] = 1

            #find min cost of path ending at N[i][j] using edge H[i][j]
            Q[i][j] = max(Q[i][j-1], R[i][j-1] + u)
            if Q[i][j] == Q[i][j-1]:
                f[i][j-1] = 1
            if Q[i][j] == R[i][j-1] + u:
                g[i][j-1] = 1

            R[i][j] = max([P[i][j], Q[i][j], R[i-1][j-1] + score(s[i-1], t[j-1])])
            if R[i][j] == P[i][j]:
                a[i][j] = 1
            if R[i][j] == Q[i][j]:
                b[i][j] = 1
            if R[i][j] == R[i-1][j-1] + score(s[i-1], t[j-1]):
                c[i][j] = 1

    #----------edge assignment-------------
    for i in reversed(range(I)):
        for j in reversed(range(J)):
            #if there is no optimal path passing through node N[i][j] which has cost R[i][j]
            #at node N[i][j], remove edges V[i][j], H[i][j] and D[i][j]
            if (a[i+1][j] == 0 or e[i][j] == 0) and (b[i][j+1] == 0 or g[i][j] == 0) and (c[i+1][j+1] == 0):
                a[i][j] = 0
                b[i][j] = 0
                c[i][j] = 0

            # if there exists optimal path passing through node N[i][j]
            if not (a[i+1][j] == 0 and b[i][j+1] == 0 and c[i+1][j+1] == 0):
                # if V[i+1][j] is an optimal path and requires edge V[i][j] to be in an optimal path, determine if an optimal path that uses edge V[i+1][j] must use edge V[i][j] and the converse:
                if a[i+1][j] == 1 and d[i][j] == 1:
                    d[i+1][j] = 1-e[i][j]
                    e[i][j] = 1-a[i][j]
                    a[i][j] = 1
                else:
                    d[i+1][j] = 0
                    e[i][j] = 0

                #if edge H[i][j+1] is in an optimal path and requires edge H[i][j] to be in an optimal path, determine if an optimal path that uses edge H[i][j+1] must use edge H[i][j] and the converse:
                if b[i][j+1] == 1 and f[i][j] == 1:
                    f[i][j+1] = 1-g[i][j]
                    g[i][j] = 1-b[i][j]
                    b[i][j] = 1
                else:
                    f[i][j+1] = 0
                    g[i][j] = 0

    #create_path(a, b, c, d, e, f, g, s, t)


    #print(R[I-1][J-1])

    return(R[I-1][J-1])
    
    print()

