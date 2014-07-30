#SYMPY IMPLEMENTATION

from sympy import *

# This function calc the phi matrix used in the AMA alg. 

def makePhi(q, cof, nlag, nlead, neq):

    # Fix size of Phi Matrix: (neq X neq)
    phi = zeros(neq, neq)

    # Calc the H_0 Matrix: (neq X neq)
    H_0 = cof[ :, (neq * nlag) : (neq * (nlag + 1) ) ]

    #print("COF")
    #pprint(cof)

    # Calc the H_+ Matrix: ( (nlead * neq) X neq )
    H_plus = cof[ :, neq * (nlag + 1) : neq * (nlag + nlead + 1) ]

    #print("H_plus")
    #pprint(H_plus)

    # Calc the Q_L Matrix: ( (neq * nlead) X (neq * nlag) )
    Q_L = q[ :, 0 : (neq * nlag)]

    # Calc the Q_R Matrix: ( (neq * nlead) X (neq * nlead) )
    Q_R = q[ :, (neq * nlag) : neq * (nlag + nlead) ]

    # Calc the B Matrix, B = (Q_R)^-1 * Q_L: (neq * nlead X neq * nlag)
    B = -1 * Q_R.inv() * Q_L

    # Calc the B_R Matrix: ( neq*nlead X neq )
    B_R = B[ :, neq * (nlag - 1) : neq * nlag]

    #print("B_R")
    #pprint(B_R)

    # Calc the phi matrix, phi = (H_0 + H_plus * B_R)^-1: (neq X neq)
    temp2 = H_0 + (H_plus * B_R)

    #print("TEMP 2")
    #pprint(temp2)


    phi = temp2.inv()

    #print("PHI")
    #pprint(phi)

    return phi
