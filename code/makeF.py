# SYMPY IMPLEMENTATION
# This function calculates the F matrix used in the AMA algorithm.

# Import the sympy package
from sympy import *

def makeF(phi, cof, q, nlag, nlead, neq):

    # Fix the size of the F Matrix: (nlead * neq X nlead * neq)
    F = zeros(nlead*neq)
    
    # Calc the H_+ Matrix: (nlead * neq X neq)
    H_plus = cof[ :, neq * (nlag + 1): neq * (nlag + nlead + 1)]

    # Calc the Q_L Matrix: (neq*nlead X neq*nlag):
    Q_L = q[:, 0 : neq*nlag]

    # Calc the Q_R Matrix: (neq*nlead X neq*nlead)
    Q_R = q[ : , neq*nlag : neq*(nlag + nlead)]

    # Calc the B Matrix, B = (Q_R)^-1 * Q_L: (neq*nlead X neq*nlag)
    B = -1 * Q_R.inv() * Q_L

    # Calc the B_R Matrix: (neq * nlead X neq)
    B_R = B[ :, neq * (nlag - 1) : neq * nlag]

    # Fill in the identity matrices of F
    # Est. the vertical index for inserting ones into F

    for i in range(0, neq * (nlead - 1) ):
        # Est. horizontal index for putting ones into F
        j = i + neq
        # Insert these identity values into F
        F[i, j] = 1

    # Construct the B_R Matrix used in the last column of F
    newB_R = zeros( (neq * nlead), neq)
    
    # Insert Identity matrix into first element of newB_R matrix
    newB_R[0 : neq, 0 : neq] = eye(neq)

    # Insert B_R^theta into newB_R for theta from 1 to theta_1
    for a in range(1, nlead):
        temp = B_R[neq * (a-1) : neq * a, :]
        newB_R[neq * a : neq * (a + 1), :] = temp
    # This will construct newB_R completely 

    # Calc the final row of matrices in F
    # Do so by looping over values raning from 1 to nlead
    for k in range(nlead, 0, -1): #extra -1 param??
        #Calc matrix to be inserted in F
        newEntry = -phi * H_plus * newB_R
        F[neq * (nlead - 1) : neq*nlead, neq * (k-1) : neq*k] = newEntry #insert in F
        if k > 1:
            # Update newB_R for next insertion matrix
            # do so by shifting each element L rows down
            for alpha in range( neq * (nlead - 1) - 1, -1, -1): #reverse order
                for beta in range(0, neq):
                    newB_R[ alpha + neq, beta] = newB_R[ alpha, beta]

            #Insert zeros in cells which are now empty
            newB_R[ 0:neq, 0:neq] = zeros(neq)

    #print("F:")
    #pprint(F)
    
    return F
    
