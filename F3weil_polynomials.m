/*
 Weil polynomials downloaded from the LMFDB on 23 October 2019.
 Below is a list (called data), collecting the weight 1 L-polynomial
 attached to each isogeny class of an abelian variety.

*/
P<x> := PolynomialRing(Integers());
F3data := [*\
9*x^4 - 18*x^3 + 15*x^2 - 6*x + 1,\
9*x^4 - 15*x^3 + 12*x^2 - 5*x + 1,\
9*x^4 - 12*x^3 + 8*x^2 - 4*x + 1,\
9*x^4 - 12*x^3 + 9*x^2 - 4*x + 1,\
9*x^4 - 12*x^3 + 10*x^2 - 4*x + 1,\
9*x^4 - 9*x^3 + 5*x^2 - 3*x + 1,\
9*x^4 - 9*x^3 + 6*x^2 - 3*x + 1,\
9*x^4 - 9*x^3 + 7*x^2 - 3*x + 1,\
9*x^4 - 9*x^3 + 8*x^2 - 3*x + 1,\
9*x^4 - 6*x^3 + x^2 - 2*x + 1,\
9*x^4 - 6*x^3 + 2*x^2 - 2*x + 1,\
9*x^4 - 6*x^3 + 3*x^2 - 2*x + 1,\
9*x^4 - 6*x^3 + 4*x^2 - 2*x + 1,\
9*x^4 - 6*x^3 + 5*x^2 - 2*x + 1,\
9*x^4 - 6*x^3 + 6*x^2 - 2*x + 1,\
9*x^4 - 6*x^3 + 7*x^2 - 2*x + 1,\
9*x^4 - 3*x^3 - 2*x^2 - x + 1,\
9*x^4 - 3*x^3 - x^2 - x + 1,\
9*x^4 - 3*x^3 - x + 1,\
9*x^4 - 3*x^3 + x^2 - x + 1,\
9*x^4 - 3*x^3 + 2*x^2 - x + 1,\
9*x^4 - 3*x^3 + 3*x^2 - x + 1,\
9*x^4 - 3*x^3 + 4*x^2 - x + 1,\
9*x^4 - 3*x^3 + 5*x^2 - x + 1,\
9*x^4 - 3*x^3 + 6*x^2 - x + 1,\
9*x^4 - 6*x^2 + 1,\
9*x^4 - 5*x^2 + 1,\
9*x^4 - 4*x^2 + 1,\
9*x^4 - 3*x^2 + 1,\
9*x^4 - 2*x^2 + 1,\
9*x^4 - x^2 + 1,\
9*x^4 + 1,\
9*x^4 + x^2 + 1,\
9*x^4 + 2*x^2 + 1,\
9*x^4 + 3*x^2 + 1,\
9*x^4 + 4*x^2 + 1,\
9*x^4 + 5*x^2 + 1,\
9*x^4 + 6*x^2 + 1,\
9*x^4 + 3*x^3 - 2*x^2 + x + 1,\
9*x^4 + 3*x^3 - x^2 + x + 1,\
9*x^4 + 3*x^3 + x + 1,\
9*x^4 + 3*x^3 + x^2 + x + 1,\
9*x^4 + 3*x^3 + 2*x^2 + x + 1,\
9*x^4 + 3*x^3 + 3*x^2 + x + 1,\
9*x^4 + 3*x^3 + 4*x^2 + x + 1,\
9*x^4 + 3*x^3 + 5*x^2 + x + 1,\
9*x^4 + 3*x^3 + 6*x^2 + x + 1,\
9*x^4 + 6*x^3 + x^2 + 2*x + 1,\
9*x^4 + 6*x^3 + 2*x^2 + 2*x + 1,\
9*x^4 + 6*x^3 + 3*x^2 + 2*x + 1,\
9*x^4 + 6*x^3 + 4*x^2 + 2*x + 1,\
9*x^4 + 6*x^3 + 5*x^2 + 2*x + 1,\
9*x^4 + 6*x^3 + 6*x^2 + 2*x + 1,\
9*x^4 + 6*x^3 + 7*x^2 + 2*x + 1,\
9*x^4 + 9*x^3 + 5*x^2 + 3*x + 1,\
9*x^4 + 9*x^3 + 6*x^2 + 3*x + 1,\
9*x^4 + 9*x^3 + 7*x^2 + 3*x + 1,\
9*x^4 + 9*x^3 + 8*x^2 + 3*x + 1,\
9*x^4 + 12*x^3 + 8*x^2 + 4*x + 1,\
9*x^4 + 12*x^3 + 9*x^2 + 4*x + 1,\
9*x^4 + 12*x^3 + 10*x^2 + 4*x + 1,\
9*x^4 + 15*x^3 + 12*x^2 + 5*x + 1,\
9*x^4 + 18*x^3 + 15*x^2 + 6*x + 1*]
;
