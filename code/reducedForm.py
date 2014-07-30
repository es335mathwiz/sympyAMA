# SYMPY IMPLEMENTATION

from sympy import *

def reducedForm(q, qrows, qcols, bcols, neq, condn):

    # Compute the reduced-form coefficient matrix, b.

    #  qs = SparseMatrix(q) ##THIS MIGHT HAVE TO CHANGE TO A LINKEDLIST FORMAT
    left = list( range( 0, qcols - qrows) )
    right = list( range( qcols - qrows, qcols) )
    
    #NONSING MIGHT NEED TO BE CHANGED, NON-EQUIVALENT SYMPY FUNCTIONS
    # nonsing = 1 / SparseMatrix(qs[:, qcols - qrows : qcols]).condition_number() > condn
    nonsing = True #assume nonsingular, invertible for symb representation
    Q_L = q[:, 0 : (qcols-qrows) ]
    Q_R = q[:, (qcols - qrows) : qcols]
    try:
        b = -1 * Q_R.inv() * Q_L
    except:
        print("Matrix is uninvertible.")

    #print("nonsing:")
    #print(nonsing)

    #print("B:")
    #pprint(b)

    return nonsing, b
