#####################################################
#E is the value of n^{-1}E(vgap) that we are looking for


import random
from GGP import Find_score

from cmath import inf
import numpy as np
import matplotlib.pyplot as plt

def randomProteins(n):
    S = ''
    T = ''
    for i in range(n):
        x = random.randint(0, 1)
        y = random.randint(0, 1)
        if x == 1:
            S += 'a'
        else:
            S += 'b'

        if y == 1:
            T += 'a'
        else:
            T += 'b'

    return S, T

#finds the E value for inputted protein length n
def E_func(n):
    u = -3
    #n = 100
    num_of_iterations = 10
    score = 0
    Total = 0
    for i in range(num_of_iterations):
        S, T = randomProteins(n)
        score = Find_score(S, T, u)
        Total += score

    print(score)
    E = Total/(num_of_iterations*n)
    print(E)
    return E


def main():
    x = [(10+20*i) for i in range(100)]
    y = [0 for i in range(100)]
    for i in range(100):
        y[i] = E_func(x[i])

    # plt.plot(x, y, 'r')
    # plt.show()

    #estimation varies more for shorter proteins, so consider more samples for shorter proteins and less samples for longer proteins?
    #Maybe so that length*(sample size) are equal for all the different lengths of proteins.
    #how do I justify this

    print()


if __name__ == '__main__':
    main()



