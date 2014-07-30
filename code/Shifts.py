#SYMPY IMPLEMENTATION
"""
Author: Alex Ahmed
Based off of files from Gary Anderson
"""
#Import sympy packages
from sympy import *
from rowSum import rowSum #own func made to sum rows of a matrix
from diagonal import diagonal

def shiftRight(x, n): #shift the rows of x to the right by n columns
                           #will leave zeros in the first n columns
                           #matrix is matrix to be augmented
                           #n is amount to shift by
    rows, cols = x.shape
    
    y = zeros(rows, cols)
    y[0:rows, n:cols] = x.extract( list(range(0,rows)), list(range(0, cols-n)))
    return y

   # V = zeros(rows, 1)
    #while( n > 0 ):
     #   m.col_insert(0,V)
      #  n = n - 1
    #return m  """


def exactShift(h,q,iq,qrows,qcols,neq):

    # h       Structural coefficient matrix ( neq, neq * ( nlag + nlead + 1 )
    # q       Q matrix
    # iq      ???
    # qrows   Number of rows in Q, neq * nlead
    # qcols   Number of columns in Q, neq * (nlag + nlead)
    # neq     Number of equations, 

    #COMPUTE THE EXACT SHIFTRIGHTS and STORE IN q. 
    
    hs = h
    nexact = 0
    left = list(range(0, qcols))
    right = list(range(qcols, qcols + neq))
    zerorows = list()
    sumVector = rowSum(hs[:, qcols: qcols+neq]) #abs built into rowSum
    sumVectorRows, sumVectorCols = sumVector.shape
    for i in range(0, sumVectorRows):
        if sumVector[i,0] == 0:
            zerorows.append(i)

    while len(zerorows) > 0 and iq <= qrows:
        nz = len(zerorows)
        hsRows, hsCols = hs.shape
        q[iq:iq+nz,0:qcols] = hs.extract(zerorows,left)
        #pprint(hs)
        #pprint(zerorows)
        #pprint(neq)
        #pprint(q)
        # pprint(hs[zerorows,:])
        #temp = hs.extract(zerorows, list(range(0, qcols+neq)))
        # hs[zerorows,:] = shiftRight( temp, neq )
        for row in zerorows:
            hs[row, :] = shiftRight( hs[row,:], neq)
        iq = iq + nz
        nexact = nexact + nz
        while len(zerorows) > 0:
            zerorows.pop()
        newSumVector = rowSum( hs[:, qcols: qcols + neq] )
        newSumVectorRows, newSumVectorCols = newSumVector.shape
        for i in range(0, newSumVectorRows):
            if newSumVector[i,0] == 0:
                zerorows.append(i)
    
    h = hs

    return h, q, iq, nexact

#########################################################################

def qrShift(h,q,iq,qrows,qcols,neq,condn):

    # Compute the numeric shiftrights and store them in q.
    
    #print(h[:, qcols: qcols+neq]) #########################

    nnumeric = 0
    left = list(range(0,qcols))
    right = list(range(qcols,qcols+neq))
    print("h:=")
    pprint(h[:,qcols:qcols+neq])
    print(h[:,qcols:qcols+neq])
    Q, R = (h[:, qcols: qcols+neq]).QRdecomposition()
    zerorows = list()
    testVector = diagonal(R)
    for i in range(0,len(testVector)):
        if testVector[i] <= condn:
            zerorows.append(i)
    
    while len(zerorows) > 0 and iq <= qrows:
        h = Q.T * h
        nz = len(zerorows)
        q[ iq: iq+nz, 0: qcols ] = h.extract(zerorows, left)
        for row in zerorows:
            hs[row, :] = shiftRight( hs[row,:], neq)
        iq = iq + nz
        nnumeric = nnumeric + nz
        Q, R  = (h[:,qcols: qcols+neq]).QRdecomposition()
        zerorows = list()
        testVector = diagonal(R)
        for i in range(0,len(testVector)):
            if testVector[i] <= condn:
                zerorows.append(i)
                
    return h, q, iq, nnumeric
