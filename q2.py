def edit_distance(s, t):
    I = len(s) + 1
    J = len(t) + 1
    D = [[0 for i in range(J)] for j in range(I)]
    
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

    return D[len(s)][len(t)]

    

def compare(a, b):
    if a == b:
        return 0
    else:
        return 1

def main():

    ed = edit_distance('shesells', 'seashells')
    print(ed)

    print()

if __name__ == '__main__':
    main()


