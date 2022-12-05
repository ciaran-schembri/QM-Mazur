//create igusa clebsch invariants with certain properties like 2-torsion.

AttachSpec("spec");

Rx<x>:=PolynomialRing(Rationals());
//this is phi : X^*(6,P2) -> X^*(6,1)
X := Curve(ProjectiveSpace(PolynomialRing(Rationals(), 2)));
KX<x> := FunctionField(X);
//phi := KX!(32/27*x^3 - 8/9*x - 8/27);

phi := KX!(-108*x^20 + 432*x^16 - 648*x^12 + 432*x^8 - 108*x^4)/(x^24 - 66*x^20 + 1023*x^16 + 2180*x^12 + 1023*x^8 - 66*x^4 + 1);
phi_perms:=[ LFT(phi,perm) : perm in Setseq(Permutations({1,2,3})) ];


height_init:=1;
small_ht:=Setseq(Set([ a/b : a,b in [-height_init..height_init] | a*b ne 0 ]));

torsion2_jsQ:=[];
primary_torsion:=[];
size:=10000;
height:=height_init;
while #torsion2_jsQ lt size do
  height:=height+1;
  extra_a:=[ height/b : b in [-height..height] | b ne 0 and GCD(height,b) eq 1 ];
  extra_b:=[ a/height : a in [-height..height] | a ne 0 and GCD(height,a) eq 1 ];
  extra := extra_a cat extra_b;
  for q in extra do
    //for phi_perm in phi_perms do
      j:=Evaluate(phi_perms[4],Place(X![q,1]));
      //Index(phi_perms,phi_perm);
      if j notin [-16/27, 0, 81/64] and Discriminant(QuaternionAlgebra< Rationals() | -6*j, -2*(27*j+16) >) eq 1 and j notin torsion2_jsQ then
        Append(~torsion2_jsQ,j);
        j; IgusaClebsch:=BabaGranathIgusaClebschDisc6(j);
        tors_heur:=TorsionGroupHeuristicUpToTwist(IgusaClebsch : bound:=100);
        invr_heur:=< PrimaryAbelianInvariants(gp) : gp in tors_heur >;
        inv:=Sprint(invr_heur);
        Append(~primary_torsion,inv);
        data:=Sprintf("%o|%o",IgusaClebsch, inv); data;
        //PrintFile("Data/X*(6,1)-curves.m",data);
      end if;
    //end for;
  end for;
end while;


/*for j in torsion2_jsQ do

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
end for;*/
