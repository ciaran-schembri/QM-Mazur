
Rx<x>:=PolynomialRing(Rationals());
//this is phi : X^*(6,P2) -> X^*(6,1)
phi := Rx!(32/27*x^3 - 8/9*x - 8/27);

height_init:=0;
small_ht:=Setseq(Set([ a/b : a,b in [-height_init..height_init] | a*b ne 0 ]));

torsion2_jsQ:=[];
size:=500;
height:=height_init;
while #torsion2_jsQ lt size do
  height:=height+1;
  extra:=[ a/b : a,b in [-height..height] | a*b ne 0 and Abs(a) eq height or Abs(b) eq height ];
  for q in extra do
    if q notin small_ht then
      Append(~small_ht,q);
      j:=Evaluate(phi,q);
      if j notin [-16/27, 0, 81/64] and Discriminant(QuaternionAlgebra< Rationals() | -6*j, -2*(27*j+16) >) eq 1 then
        Append(~torsion2_jsQ,j);
      end if;
    end if;
  end for;
end while;


for j in torsion2_jsQ do
  IgusaClebsch:=BabaGranathIgusaClebschDisc6(j);
  tors2_heur:=TorsionHeuristicUpToTwistDivisibleBy(IgusaClebsch);
  Write("Data/X*(6,P2)-curves.m", Sprintf("%o|%o",IgusaClebsch, tors2_heur));
end for;

end for;
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
