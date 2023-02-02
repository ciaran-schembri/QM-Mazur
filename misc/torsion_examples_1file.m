AttachSpec("../CHIMP/CHIMP.spec");
Rx<x>:=PolynomialRing(Rationals());
prec := 500;
F := RationalsExtra(prec);
///////////////////////////
//(Z/2Z)^2
f:=-180*x^6 - 159*x^5 + 894*x^4 + 1691*x^3 + 246*x^2 - 672*x + 80;
C:=HyperellipticCurve(f);
J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [2,2];
//T; (Z/2)^2;

X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : Geometric:=true);
assert IsQuaternionAlgebra(B);

_,E:=HeuristicEndomorphismAlgebra(X : Geometric:=false);
assert IsCommutative(E);
b:=Basis(E)[2];
MinimalPolynomial(b);
assert Discriminant(NumberField(MinimalPolynomial(b))) eq 12;



////////////////////////////////
//(Z/3)^2
C:=HyperellipticCurve([Polynomial([RationalField() | 105, 270, -45, -270, 315, -270, -15]), Polynomial([RationalField() |])]);
J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [3,3];
//T; (Z/3)^2;

X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : Geometric:=true);
assert IsQuaternionAlgebra(B);

_,E:=HeuristicEndomorphismAlgebra(X : Geometric:=false);
assert Dimension(E) eq 2;
assert IsCommutative(E);
b:=Basis(E)[2];
MinimalPolynomial(b);
assert Discriminant(NumberField(MinimalPolynomial(b))) eq 8;

//R:=HeuristicEndomorphismRepresentation(X : Geometric:=false);


/////////////////////////////////////
//(Z/3)
C:=HyperellipticCurve([Polynomial([RationalField() | -560, -1200, 120, 800, 420,
-300, -10]), Polynomial([RationalField() |])]);

J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [3];
//T; (Z/3);

X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : Geometric:=true);
assert IsQuaternionAlgebra(B);

_,E:=HeuristicEndomorphismAlgebra(X : Geometric:=false);
assert Dimension(E) eq 1;



////////////////////////////////////
//Z/6Z
C:=HyperellipticCurve([Polynomial([RationalField() | -343, 0, 294, -49, -63, 21, 5]), Polynomial([RationalField() |])]);

J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [2,3];
//T; (Z/3);

X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : Geometric:=true);
assert IsQuaternionAlgebra(B);

_,E:=HeuristicEndomorphismAlgebra(X : Geometric:=false);
assert Dimension(E) eq 1;


//////////////////////////////////////
//2-torsion discriminant 15 example
C:=HyperellipticCurve([Polynomial([RationalField() | 25039, 74970, 47271, -14780, 23673, 3930, -367]), Polynomial([RationalField() |])]);
J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [2];


prec := 1000;
F := RationalsExtra(prec);
X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : CC:=true);
B;
_,D:= IsQuaternionAlgebra(B);
Discriminant(D); //15

/////////////////////////////////////////

IgusaClebsch:=ChangeUniverse([ 36192, -567000000, -4136184000000, 93312000000000000 ],Rationals());
X:=HyperellipticCurveFromIgusaClebsch(IgusaClebsch);
C:=ReducedWamelenModel(X);
tors_heur:=TorsionGroupHeuristicUpToTwist(IgusaClebsch : bound:=400);
J:=Jacobian(C);
T:=TorsionSubgroup(J);
