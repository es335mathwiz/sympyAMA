#SYMPY IMPLEMENTATION
from sympy import *

def colSum(matrix):
    sumVector = zeros(1, matrix.cols)
    for col in range(0, matrix.cols):
        for row in range(0, matrix.rows):
            sumVector[0, col] = sumVector[0, col] + abs(matrix[row, col])
    return sumVector
