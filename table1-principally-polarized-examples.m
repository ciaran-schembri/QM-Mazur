//Need to have https://github.com/edgarcosta/endomorphisms to independently check the endomorphisms algebras,
//however, the geometric endomorphism algebra is already known because of the Igusa-Invariants in Lin-Yang (and Baba-Granath).

Rx<x>:=PolynomialRing(Rationals());
prec := 1000;
F := RationalsExtra(prec);

//////////////////////////////
//line1: A(Q)_tors=(Z/2Z) and D =10
C:=HyperellipticCurve([Polynomial([RationalField() | 0, -1312695, 0, 2187825, 0, -729275, -145855]), Polynomial([RationalField() |])]);
J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [2];

X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : Geometric:=true);
assert IsQuaternionAlgebra(B);
tr,D:=IsQuaternionAlgebra(B);
assert Discriminant(D) eq 10;


_,E:=HeuristicEndomorphismAlgebra(X : Geometric:=false);
assert Dimension(E) eq 1;


//////////////////////////////
//line2: A(Q)_tors = (Z/2Z)^2 and D=6
f:=-180*x^6 - 159*x^5 + 894*x^4 + 1691*x^3 + 246*x^2 - 672*x + 80;
C:=HyperellipticCurve(f);
J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [2,2];

X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : Geometric:=true);
assert IsQuaternionAlgebra(B);
tr,D:=IsQuaternionAlgebra(B);
assert Discriminant(D) eq 6;


_,E:=HeuristicEndomorphismAlgebra(X : Geometric:=false);
assert IsCommutative(E);
b:=Basis(E)[2];
MinimalPolynomial(b);
assert Discriminant(NumberField(MinimalPolynomial(b))) eq 12;


//////////////////////////////////////////
//line3: A(Q)_tors = Z/3Z and  D=15
C:=HyperellipticCurve([Polynomial([RationalField() | -634465, -540930, -43680, 234260, 602160, 345930, 17095]), Polynomial([RationalField() |])]);
J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [3];

X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : Geometric:=true);
assert IsQuaternionAlgebra(B);
tr,D:=IsQuaternionAlgebra(B);
assert Discriminant(D) eq 15;

_,E:=HeuristicEndomorphismAlgebra(X : Geometric:=false);
assert Dimension(E) eq 1;



////////////////////////////////
//line4: A(Q)_tors = (Z/3Z)^2 and D=6
C:=HyperellipticCurve([Polynomial([RationalField() | 105, 270, -45, -270, 315, -270, -15]), Polynomial([RationalField() |])]);
J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [3,3];

X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : CC:=true);
assert IsQuaternionAlgebra(B);
tr,D:=IsQuaternionAlgebra(B);
assert Discriminant(D) eq 6;


_,E:=HeuristicEndomorphismAlgebra(X : Geometric:=false);
assert Dimension(E) eq 2;
assert IsCommutative(E);
b:=Basis(E)[2];
MinimalPolynomial(b);
assert Discriminant(NumberField(MinimalPolynomial(b))) eq 8;


////////////////////////////////////
//line5: A(Q)_tors = (Z/6Z) and D=6Z
C:=HyperellipticCurve([Polynomial([RationalField() | -343, 0, 294, -49, -63, 21, 5]), Polynomial([RationalField() |])]);
J:=Jacobian(C);
T:=TorsionSubgroup(J);
assert PrimaryAbelianInvariants(T) eq [2,3];

X := ChangeRing(C,F);
_,B:=HeuristicEndomorphismAlgebra(X : Geometric:=true);
assert IsQuaternionAlgebra(B);
tr,D:=IsQuaternionAlgebra(B);
assert Discriminant(D) eq 6;

_,E:=HeuristicEndomorphismAlgebra(X : Geometric:=false);
assert Dimension(E) eq 1;

