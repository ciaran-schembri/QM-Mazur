//////////////////////////////////////////////

load "QM-Mazur/misc/prop4-3-9-Z2Z2-defining-equations.m";


Rx<x>:=PolynomialRing(Rationals());
//this is phi : X^*(6,P2) -> X^*(6,1)
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
C:=HyperellipticCurveFromIgusaInvariants(BG6Igusa);









