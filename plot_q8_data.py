with open('q8_data.txt', 'r') as f:
    data = [[num for num in line.split()] for line in f]

import matplotlib.pyplot as plt
import numpy as np

x = [0 for i in range(len(data))]
y = [0 for i in range(len(data))]

for i in range(len(data)):
    x[i] = float(data[i][0])
    y[i] = float(data[i][1])

plt.scatter(x, y, color = 'r', s = 1)
#plt.plot(x, y)

plt.axhline(y = 0.433, color = 'b', linewidth = 0.8, linestyle = 'dashed')

# plt.plot(x, np.poly1d(np.polyfit(x, y, 10))(x))

plt.title('Graph of n against estimated E')
plt.xlabel('n')
plt.ylabel('E')

plt.savefig('q8_1.jpg', dpi=300)
plt.show()



