#SYMPY IMPLEMENTATION

#REDO/newer

#import sympy packages
from sympy import *
from colSum import colSum

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

    sumVector = colSum(a)
    sumVectorRows, sumVectorCols = sumVector.shape
    for i in range(0, sumVectorCols):
        if sumVector[0, i] == 0:
            zerocols.append(i)

    while len(zerocols) > 0:

        ##
        zerocols.reverse() ##reverse order to preserve indicies 
        for n in zerocols:
            a.col_del(n)   #eliminate zero columns
        for n in zerocols:
            a.row_del(n)   #eliminate zero rows
        for z in zerocols:
            del js[z]      #update list js
        zerocols.reverse() #put back in order
        while len(zerocols) > 0:
            zerocols.pop()
        sumVector2 = colSum(a)
        sumVector2Rows, sumVector2Cols = sumVector2.shape
        for i in range(0, sumVector2Cols):
            if sumVector2[0,i] == 0:
                zerocols.append(i)
    ia = len(js)
    return a, ia, js

   
