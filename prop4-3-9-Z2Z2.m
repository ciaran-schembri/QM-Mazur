//////////////////////////////////////////////
//this loads the hyperelliptic polynomial which defines the genus 2 curves with Igusa invariants below.
load "QM-Mazur/misc/prop4-3-9-Z2Z2-defining-equations.m";

X := Curve(ProjectiveSpace(PolynomialRing(Rationals(), 2)));
KX<x> := FunctionField(X);

j:=(-64*x^20 + 256*x^16 - 384*x^12 + 256*x^8 - 64*x^4)/(x^24 + 42*x^20 + 591*x^16 + 2828*x^12 + 591*x^8 + 42*x^4 + 1);
//j:=phi;
J2:=12*(j+1);
J4:=6*(j^2+j+1);
J6:=4*(j^3-2*j^2+1);
J8:=(J2*J6-J4^2)/4;
J10:=j^3;
BG6Igusa:=[J2,J4,J6,J8,J10];

//this shows that End(A_Q) = Q(sqrt(3))
assert ((x^12 - 33*x^8 - 33*x^4 + 1)/(2*x^10 - 4*x^6 + 2*x^2))^2 eq -27 - 16/j;

//this shows that the mestre obstruction is generically split.
a1:=6/(x^8 + 14*x^4 + 1);
b1:=(8*x^10 - 16*x^6 + 8*x^2)/(x^8 + 14*x^4 + 1);
assert a1*b1^2 eq -6*j;

a2:=-2/(x^8 + 14*x^4 + 1);
b2:=(-(4*x^12 - 132*x^8 - 132*x^4 + 4)/(x^8 + 14*x^4 + 1));
assert a2*b2^2 eq -2*(27*j+16);

assert Discriminant(QuaternionAlgebra< Rationals() | -2,6 >) eq 1;










