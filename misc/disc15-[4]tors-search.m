


Rx<x>:=PolynomialRing(Rationals());
//this is phi : X^*(6,P2) -> X^*(6,1)
X := Curve(ProjectiveSpace(PolynomialRing(Rationals(), 2)));
KX<x> := FunctionField(X);
phi1Elkies := KX!(-1*(16*x^6-24*x^4+9*x^2-1))-1;
phi1LY:=-27*phi1Elkies;

//X^*(15,4) -> X^*(15,1) belyi map

