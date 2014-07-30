#SYMPY IMPLEMENTATION
#import sympy packages
from sympy import *

def copyW(q, w, js, iq, qrows):

    #Copy eigenvectors corresponding to largest roots onto remaining
    #empty rows and columns (js) of q

    if iq < qrows:
        lastrows = list(range(iq, qrows))
        wrows = list(range(0, len(lastrows)))

        #print("SPLICED W")
        #pprint(w[:, 0 : len(lastrows)].T)

        #print("q")
        #pprint(q)

        #print("q splice")
        #pprint(q[iq : qrows, 2 : 3 ])
        jsend = js[len(js) - 1]
        q[iq : qrows+1, js[0] : jsend+1 ] = w[:, 0 : len(lastrows)].T
        #for c in js:
            
        #print("THIS IS JS:")
        #pprint(js)

        #print("FINAL Q")
        #pprint(q)
        
        #check list boundaries for splicing, +1 may be necessary 
        #return spliced submatrix

    return q
