import numpy as np
csv_name = "matrix.csv"

r = np.genfromtxt(csv_name, delimiter=',', names=True, case_sensitive=True)
print(repr(r))
for i in range(len(r)):
    for j in range(len(r)):
        print(str(abs(r[i][j] - r[j][i]))+" "+str(i)+' '+str(j))

