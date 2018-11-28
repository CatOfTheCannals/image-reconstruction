import numpy as np

def load_cpp_matrix(mat_name):
    mat_path = "./debug_matrices/{}.mat".format(mat_name)
    file = open(mat_path, "r") 
    rows = int(file.readline())
    cols = int(file.readline())
    mat = np.empty((rows,cols))
    k = 0
    while True:
        number = file.readline()
        if number == '':
            break
        i = k // cols
        j = k % cols
        mat[i][j] = float(number)    
        k += 1
    return mat

def save_matrix_to_cpp(matrix, mat_name):
    mat_path = "./debug_matrices/{}.mat".format(mat_name)
    file = open(mat_path, "w") 
    if len(matrix.shape) == 1:
        matrix = matrix.reshape(-1,1)
    rows = matrix.shape[0]
    cols = matrix.shape[1]
    file.write('{}\n'.format(rows))
    file.write('{}\n'.format(cols))

    for i in range(rows):
        for j in range(cols):
            file.write('{}\n'.format(matrix[i][j]))

