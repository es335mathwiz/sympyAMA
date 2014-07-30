#SYMPY IMPLEMENTATION
from sympy import *

from Shifts import exactShift, qrShift
from buildA import buildA
from Eigensystem import Eigensystem
from copyW import copyW
from makePhi import makePhi
from makeF import makeF
from reducedForm import reducedForm

#this file is intended for testing purposes
#used to initialize values necessary to test

def initValues():

    #define symbols to be used
    r,d = symbols('r d')
    
    h = Matrix([[0,0,1 + r,0,-1,-1],[0,1 - d,0,1,0,0]])
    neq = 2
    nlag = 1
    nlead = 1
    condn = 0
    uprbnd = 1 
    nexact = 0
    nnumeric = 0
    lgroots = 0
    iq = 0
    aimcode = 0
    b = 0
    qrows = neq*nlead
    qcols = neq*(nlag+nlead)
    bcols = neq*nlag
    q = zeros(qrows,qcols)
    rts = zeros(qcols,1)
    #originalH = h
    originalH = Matrix([[0,0,1 + r,0,-1,-1],[0,1 - d,0,1,0,0]])  #??? why necessary

    print("OG_H initial")
    pprint(originalH)
    
    h, q, iq, nexact = exactShift(h, q, iq, qrows, qcols, neq)

    print("H after exactShift")
    pprint(h)
    print("Q after exactShift")
    pprint(q)

    h, q, iq, nnumeric = qrShift(h, q, iq, qrows, qcols, neq, condn)

    print("H after qrShift")
    pprint(h)
    print("Q after qrShift")
    pprint(q)

    a, ia, js = buildA(h, qcols, neq)

    print("A after buildA")
    pprint(a)

    print("ia:")
    pprint(ia)

    print("js:")
    pprint(js)

    w, rts = Eigensystem(a, uprbnd, min(len(js), qrows - iq + 1))

    print("w after EIGENSYSTEM")
    pprint(w)

    print("rts:")
    pprint(rts)

    q = copyW(q, w, js, iq, qrows)

    #print("OG_H after copyW")
    #pprint(originalH)

    phi = makePhi(q, originalH, nlag, nlead, neq)
    F = makeF(phi, originalH, q, nlag, nlead, neq)
    nonsing, b = reducedForm(q, qrows, qcols, bcols, neq, condn)

    return nonsing, b
initValues()
     
