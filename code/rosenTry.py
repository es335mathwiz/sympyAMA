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

def t1():

    #define symbols to be used
    a = Symbol('a', real = True, positive = True)
    d = Symbol('d', real = True, positive = True)
    b = Symbol('b', real = True, positive = True)
    delta = Symbol('delta', real = True, positive = True)
    
    h = Matrix([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -a, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, -1 + delta, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, d, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
  0, -b *(1 - delta), 0, 0, 0, 0, 0, 0, 0, 0, 
  0, (-b**4)*(1 - delta)**4, 0, 0, 0]]) 
    neq = 4
    nlag = 4
    nlead = 4
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
    originalH = Matrix([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -a, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, -1 + delta, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, d, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
  0, -b *(1 - delta), 0, 0, 0, 0, 0, 0, 0, 0, 
  0, (-b**4)*(1 - delta)**4, 0, 0, 0]])

    #print("OG_H initial")
    #pprint(originalH)
    
    h, q, iq, nexact = exactShift(h, q, iq, qrows, qcols, neq)

    #print("OG_H after exactShift")
    #pprint(originalH)

    h, q, iq, nnumeric = qrShift(h, q, iq, qrows, qcols, neq, condn)

    #print("OG_H after qrShift")
    #pprint(originalH)

    a, ia, js = buildA(h, qcols, neq)

    #print("OG_H after buildA")
    #pprint(originalH)

    w, rts = Eigensystem(a, uprbnd, min(len(js), qrows - iq + 1))

    #print("OG_H after EIGENSYSTEM")
    #pprint(originalH)

    q = copyW(q, w, js, iq, qrows)

    #print("OG_H after copyW")
    #pprint(originalH)

    phi = makePhi(q, originalH, nlag, nlead, neq)
    F = makeF(phi, originalH, q, nlag, nlead, neq)
    nonsing, b = reducedForm(q, qrows, qcols, bcols, neq, condn)

    pprint(b)

    return nonsing, b

t1()

