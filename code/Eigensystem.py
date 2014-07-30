#SYMPY IMPLEMENTATION
from sympy import *

def Eigensystem(a, uprbnd, rowsLeft):

    #  Compute the roots and the left eigenvectors of the companion
    #  matrix, sort the roots from large-to-small, and sort the
    #  eigenvectors conformably.  Map the eigenvectors into the real
    #  domain. Count the roots bigger than uprbnd.
    
    #print("A:")
    #pprint(a)

    #print("ROWSLEFT:")
    #print(rowsLeft)


    #use transpose of a
    
    temprts = a.T.eigenvals().keys()

    ##
    #pprint(rts)
    print("a transpose:")
    pprint(a.T)

    M = a.T.eigenvects() #M is temp holder of eigenvalues in raw form


    ##
    #print("This is M")
    #pprint(M)
    temp = [ x[2] for x in M] #must eigenvectors be normalized?
    

    temp.reverse() ###JUST CHANGED
    #print("THIS IS TEMP")
    #pprint(temp)
    w = temp[0][0]
    #print("THIS IS W:")
    #pprint(w)
    for mat in temp:
        #print("THIS IS MAT:")
        #pprint(mat)
        for q in mat:

            #print("THIS IS Q:")
            #pprint(q)
            w = w.row_join(q)
            
            #print("W AFTER ROW JOIN")
            #pprint(w)
    w.col_del(0)
    #print("FINAL W")
    #pprint(w)
      

    #count = len(temp)
    #w = temp[0]
    #while count <= len(temp):
        #w = w.row_join(temp[count])
        #count = count + 1
    
    #LGROOTS MUST STILL BE DEFINED ###
    #LOOK AT MATHEMATICA FILE FOR PARAMETERS TO SUB
    rts = list()
    for val in temprts:
        rts.append(val)

    return w, rts    #eliminated lgroots because unapplicable 



