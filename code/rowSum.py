#SYMPY IMPLEMENTATION
from sympy import *

def rowSum(matrix):
    sumVector = zeros(matrix.rows, 1)
    for row in range(0,matrix.rows):
        for col in range(0, matrix.cols):
            sumVector[row,0] = sumVector[row,0] + abs(matrix[row,col])

    return sumVector
    
