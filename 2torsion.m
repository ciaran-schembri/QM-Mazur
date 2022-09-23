
Rx<x>:=PolynomialRing(Rationals());
phi := Rx!(32/27*x^3 - 8/9*x - 8/27);

height:=100;
small_ht:=Setseq(Set([ a/b : a,b in [-height..height] | a*b ne 0 ]));
torsion2_js_init:= Setseq(Set([ Evaluate(phi,a) : a in small_ht ]));
//torsion2_js_goodred:= [ a : a in torsion2_js_init | a ne 0 and (Set(PrimeDivisors(Numerator(a)*Denominator(a))) subset {2,3} eq true) ];
torsion2_jsQ := [ j : j in torsion2_js_init | not(j in [-16/27, 0, 81/64]) and Discriminant(QuaternionAlgebra< Rationals() | -6*j, -2*(27*j+16) >) eq 1 ];

for j in torsion2_jsQ do
  IgusaClebsch:=BabaGranathIgusaClebschDisc6(j);
  tors2_heur:=TorsionHeuristicUpToTwistDivisibleBy(2,IgusaClebsch : bound:=150);
  if tors2_heur notin [1,2] then
    j; tors2_heur;
  end if;
end for;


C:=HyperellipticCurveFromIgusaClebsch(IgusaClebsch);
TwoTorsionSubgroup(Jacobian(C));

for d in [-100000..100000] do
  if d ne 0 then
    if IsSquarefree(d) and d ne 0 then
      Xd:=QuadraticTwist(X,d);
      T:=TwoTorsionSubgroup(Jacobian(Xd));
      if #T ne 2 then
        d; T;
      end if;
    end if;
  end if;
end for;
