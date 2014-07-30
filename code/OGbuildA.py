#SYMPY IMPLEMENTATION

#import sympy packages
from sympy import *

def buildA(h, qcols, neq):

    # Build the companion matrix, deleting inessential lags.
    # Solve for x_{t+nlead} in terms of x_{t+nlag},...,x_{t+nlead-1}.

    #Define variables
    left = list(range(0, qcols))
    right = list(range(qcols, qcols + neq)) 
    hs = SparseMatrix(h)
    hs_left = hs[:,0:qcols]
    hs_right = hs[:,qcols: qcols + neq]
    a0 = hs[:, qcols: (qcols + neq)] #might need last item in range to be + 1
    hs[:, 0: qcols] = -1 * hs_right.LUsolve(hs_left) #LUsolve left hand side of hs

    #Build big transition matrix

    a = zeros(qcols, qcols)

    if qcols > neq:
        a[0:qcols-neq, neq:qcols] = eye(qcols-neq)

    a[ (qcols - neq) : qcols, :] = hs[:, 0:qcols]

    #  Delete inessential lags and build index array js.  js indexes the
    #  columns in the big transition matrix that correspond to the
    #  essential lags in the model.  They are the columns of q that will
    #  get the unstable left eigenvectors. 

    js = list(range(0, qcols))
    zerocols = list()

    sumVector = zeros(qcols,1)

    #must use loop to sum matrix rows
    #sympy has no built in func for this
    for y in range(0, qcols):
        for x in range(0, qcols):
            sumVector[y] = sumVector[y] + abs( a[y, x] )
    #assign len
    sumVectorRows = sumVector.rows
    sumVectorCols = sumVector.cols

    for i in range(0, sumVectorCols):
        if sumVector[0, i] == 0:
            zerocols.append(i)

    #row, column annihilation 

    while len(zerocols) > 0:
        for n in zerocols:
            a.col_del(n)
            a.row_del(n)
            #delete entry from list js

            zerocols.reverse() #reverse order of zerocols to preserve indices of js
            for z in zerocols:
                # del snaq[z]
                del js[z]
            zerocols.reverse() #put back in original order
            while len(zerocols) > 0:
                zerocols.pop()

            #get updated dimensions of a

            a_rows, a_cols = a.shape

            #create sumVector2
            
            sumVector2 = zeros(a_rows,1)
            for y in range(0, a_rows):
                for x in range(0, a_cols):
                    sumVector2[y] = sumVector2[y] + abs( a[y, x] )
            sumVector2Rows, sumVector2Cols = sumVector2.shape
            for i in range(0, sumVector2Cols):
                if sumVector2[0, i] == 0:
                    zerocols.append(i)

    ia = len(js)
    return a, ia, js
