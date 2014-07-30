#SYMPY IMPLEMENTATION
from sympy import *
#takes in a matrix and returns the abs diagonal entries

def diagonal(m): #m is some matrix(usually square)
    n = m.rows
    V = zeros(n,1)
    count = 0
    while count < n:
        V[count,0] = abs( m[count,count] )
        count = count + 1
    return V
    
