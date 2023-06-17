#derived from source [1], a better version of Gotoh and Taylor's algorithm with O(mn) complexity and affine gap costs
#uses many different edge options to map out path
#The original algorithm takes w_k = v + uk, here we would take u to be 0, and v to be some fixed value
#This algorithm is called SS-2


from cmath import inf
import numpy as np

def score(x, y):
    indexx = blosum[0].index(str(x))
    indexy = blosum[0].index(str(y))
    return int(blosum[indexx][indexy])


with open('blosum.txt', 'r') as f:
    blosum = [[num for num in line.split()] for line in f]


file = open('proteins.txt', "r")
content = file.readlines()
content = [line.rstrip() for line in content]
file.close()
proteinA = content[1][0:9]
proteinB = content[3][0:9]

u = -12 #constant gap penalty
s = proteinA
t = proteinB
J = len(t)+1
I = len(s)+1

#[1]
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
        #[2]: find max cost of path ending at N[i][j] using edge V[i][j]
        P[i][j] = max([P[i-1][j], R[i-1][j] + u])
        #[3]:
        if P[i][j] == P[i-1][j]:
            d[i-1][j] = 1

        if P[i][j] == R[i-1][j] + u:
            e[i-1][j] = 1

        #[4]: find max cost of path ending at N[i][j] using edge H[i][j]
        Q[i][j] = max(Q[i][j-1], R[i][j-1] + u)
        #[5]:
        if Q[i][j] == Q[i][j-1]:
            f[i][j-1] = 1
        if Q[i][j] == R[i][j-1] + u:
            g[i][j-1] = 1

        #[6]:
        R[i][j] = max([P[i][j], Q[i][j], R[i-1][j-1] + score(s[i-1], t[j-1])])
        #[7]:
        if R[i][j] == P[i][j]:
            a[i][j] = 1
        if R[i][j] == Q[i][j]:
            b[i][j] = 1
        if R[i][j] == R[i-1][j-1] + score(s[i-1], t[j-1]):
            c[i][j] = 1

#----------edge assignment-------------
for i in reversed(range(I)):
    for j in reversed(range(J)):
        #[8]: if there is no optimal path passing through node N[i][j] which has cost R[i][j]
        #at node N[i][j], remove edges V[i][j], H[i][j] and D[i][j]
        if (a[i+1][j] == 0 or e[i][j] == 0) and (b[i][j+1] == 0 or g[i][j] == 0) and (c[i+1][j+1] == 0):
            a[i][j] = 0
            b[i][j] = 0
            c[i][j] = 0

        #[9]: if there exists optimal path passing through node N[i][j]
        if not (a[i+1][j] == 0 and b[i][j+1] == 0 and c[i+1][j+1] == 0):
            #[10]: if V[i+1][j] is an optimal path and requires edge V[i][j] to be in an optimal path, determine if an optimal path that uses edge V[i+1][j] must use edge V[i][j] and the converse:
            if a[i+1][j] == 1 and d[i][j] == 1:
                d[i+1][j] = 1-e[i][j]
                e[i][j] = 1-a[i][j]
                a[i][j] = 1
            else:
                d[i+1][j] = 0
                e[i][j] = 0

            #[11]: if edge H[i][j+1] is in an optimal path and requires edge H[i][j] to be in an optimal path, determine if an optimal path that uses edge H[i][j+1] must use edge H[i][j] and the converse:
            if b[i][j+1] == 1 and f[i][j] == 1:
                f[i][j+1] = 1-g[i][j]
                g[i][j] = 1-b[i][j]
                b[i][j] = 1
            else:
                f[i][j+1] = 0
                g[i][j] = 0

print(R[I-1][J-1])

input()