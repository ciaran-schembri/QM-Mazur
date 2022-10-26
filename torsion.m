SetColumns(0);

intrinsic PointsCardinalityTwistsFiniteField(IgusaClebschModp::SeqEnum) -> SeqEnum
  {Given Igusa-Clebsch invariants over a finite field F_p,
   create a genus 2 curve from these invariants defined over F_p.
    Return the point count of the Jacobian of this curve and it's twist}

  Rx<x>:=PolynomialRing(Integers());
  assert Characteristic(Universe(IgusaClebschModp)) gt 0;
  Cp:=HyperellipticCurveFromIgusaClebsch(IgusaClebschModp);
  assert BaseRing(Cp) eq Universe(IgusaClebschModp);
  Jp:=Jacobian(Cp);
  fp:=Rx!EulerFactor(Jp);
  coefs:=Coefficients(fp);
  if coefs[2] eq 0 and coefs[4] eq 0 then
    gp := coefs[1] - coefs[3]*x^2 + coefs[5]*x^4;
    return [ Evaluate(fp,1), Evaluate(gp,1) ];
  else
    return [ Evaluate(fp,1),Evaluate(fp,-1) ];
  end if;

end intrinsic;



intrinsic IntersectAbelianGroups(A::GrpAb,B::GrpAb) -> GrpAb
  {return the subgroup G of finite abelian groups A and B such that all other subgroups of both are contained in G}
  subgroupsA:=[ G`subgroup : G in Subgroups(A) ];
  subgroupsB:=[ G`subgroup : G in Subgroups(B) ];

  invariantsA := [ PrimaryAbelianInvariants(G) : G in subgroupsA ];
  invariantsB := [ PrimaryAbelianInvariants(G) : G in subgroupsB ];

  intersect:= Setseq(Set([ invar : invar in invariantsA | invar in invariantsB ]));
  maximal_intersects:=[];
  for H in intersect do
    intersect_exclH := Exclude(intersect,H);
    allsubs_exclH :=Setseq(Set(&cat[ [ PrimaryAbelianInvariants(C`subgroup) : C in Subgroups(AbelianGroup(D)) ] : D in intersect_exclH ]));
    if exists(t){ G : G in allsubs_exclH | H eq G } eq false then
      Append(~maximal_intersects,H); //H;
    end if;
  end for;
  //maximal_intersects;

  assert #maximal_intersects eq 1;
  return AbelianGroup(maximal_intersects[1]);
end intrinsic;

  /*end for;
  sort:=function(I,J);
    if I subset J then
      return 1;
    elif J subset I then
      return -1;
    else
      assert true eq false;
    end if;
  end function;
  Sort(intersect, sort);
  maximal_intersect:=AbelianGroup(intersect[#intersect]);*/






intrinsic JacobianGroupTwistsFiniteField(IgusaClebschModp::SeqEnum) -> SeqEnum
  {Given Igusa-Clebsch invariants over a finite field F_p,
   create a genus 2 curve from these invariants defined over F_p.
    Return the point count of the Jacobian of this curve and it's twist}

  assert Characteristic(Universe(IgusaClebschModp)) gt 0;
  Cp:=HyperellipticCurveFromIgusaClebsch(IgusaClebschModp);
  twists:=Twists(Cp);
  return [ AbelianGroup(Jacobian(X)) : X in twists ];

end intrinsic;


intrinsic BabaGranathIgusaClebschDisc6(j::FldRatElt) -> SeqEnum
  {Given j from Baba-Granath's discriminant 6 family,
  create the Igusa-Clebsch invariants of the curve associated to j.}

  // The Igusa invariants from Baba-Granath for D=6
  J2:=12*(j+1);
  J4:=6*(j^2+j+1);
  J6:=4*(j^3-2*j^2+1);
  J8:=(J2*J6-J4^2)/4;
  J10:=j^3;
  BG6Igusa:=[J2,J4,J6,J8,J10];

  //The Igusa-Clebsch invariants are
  I2:=8*J2;
  I4:=-96*J4+4*J2^2;
  I6:=-(576*J6-8*J2^3+160*J2*J4);
  I10:=4096*J10;
  BG6IC:=[I2,I4,I6,I10];

  return BG6IC;
end intrinsic;



intrinsic TorsionGroupHeuristicUpToTwist(IgusaClebsch::SeqEnum:group:=AbelianGroup([1]), bound:=150) -> RngIntElt
  {Given a Igusa-Clebsch invariants which define a curve over Q, check the biggest possible prime power torsion of the Jacobian.
  This is done by reducing mod p for p up to a bound and up to twist. }

  group_invs:=PrimaryAbelianInvariants(group);
  if Universe(IgusaClebsch) eq Integers() then
    ChangeUniverse(~IgusaClebsch,Rationals());
  end if;
  assert Universe(IgusaClebsch) eq Rationals();
  badprimes:=PrimeDivisors(&*([ Denominator(I) : I in IgusaClebsch ] cat [Numerator(IgusaClebsch[4])])*30);
  primes:=[ a : a in PrimesUpTo(bound) | a notin badprimes ];
  possible_groups:=[];
  IgusaClebschModp:=ChangeUniverse(IgusaClebsch,FiniteField(primes[1]));
  all_possible_groups:=JacobianGroupTwistsFiniteField(IgusaClebschModp);
  invs:=Setseq(Set([ PrimaryAbelianInvariants(G) : G in all_possible_groups ]));
  all_possible_groups:=[ AbelianGroup(I) : I in invs ];

  flag:=true;
  for p in primes do
    IgusaClebschModp:=ChangeUniverse(IgusaClebsch,FiniteField(p));
    group_twists:=JacobianGroupTwistsFiniteField(IgusaClebschModp);
    grps:=&cat[ [ IntersectAbelianGroups(A,B) : B in group_twists ] : A in all_possible_groups ];
    invs:=Setseq(Set([ PrimaryAbelianInvariants(G) : G in grps ]));
    all_possible_groups:=[ AbelianGroup(I) : I in invs ];
    if group_invs notin Setseq(Set(&cat[ [ PrimaryAbelianInvariants(C`subgroup) : C in Subgroups(AbelianGroup(D)) ] : D in all_possible_groups ])) then
      flag:=false;
      break p;
    end if;
  end for;

  if flag eq false then
    return "group not in torsion";
  else
    return all_possible_groups;
  end if;

end intrinsic;





intrinsic TorsionHeuristicUpToTwistDivisibleBy(IgusaClebsch::SeqEnum : bound:=200, primes:=[], divisibleby:=-1) -> RngIntElt
  {Given a Igusa-Clebsch invariants which define a curve over Q, check highest power of q that divides
  the cardinality of the point count of the Jacobian. This is done by reducing mod p
  for p up to a bound and up to twist. }

  if Universe(IgusaClebsch) eq Integers() then
    ChangeUniverse(~IgusaClebsch,Rationals());
  end if;
  assert Universe(IgusaClebsch) eq Rationals();
  badprimes:=PrimeDivisors(&*([ Denominator(I) : I in IgusaClebsch ] cat [Numerator(IgusaClebsch[4])])*30);
  all_pointcounts:=[];
  for p in [ a : a in PrimesUpTo(bound) | a notin badprimes ] do
    IgusaClebschModp:=ChangeUniverse(IgusaClebsch,FiniteField(p));
    pointcount_twists:=PointsCardinalityTwistsFiniteField(IgusaClebschModp);
    if divisibleby eq -1 or #primes ne 1 then
      Append(~all_pointcounts,pointcount_twists);
    else
      max:=primes[1]^(Integers()!Maximum([ Valuation(a,primes[1]) : a in pointcount_twists ]));
      if IsDivisibleBy(max,divisibleby) then
        Append(~all_pointcounts,pointcount_twists);
      else
        return 0;
      end if;
    end if;
  end for;

  if primes eq [] then
    primes:=PrimeDivisors(&*(all_pointcounts[1]));
  end if;
  tors:=1;
  for q in primes do
    pointcount_modp:=[];
    for point_count in all_pointcounts do
      max:=Maximum([ Valuation(a,q) : a in point_count ]);
      if max eq 0 then
        Append(~pointcount_modp,max);
        break;
      else
        //[q,max];
        Append(~pointcount_modp,max);
      end if;
    end for;
    min_valuation:=Minimum(pointcount_modp);
    tors:=tors*q^min_valuation;
  end for;
  return tors;

end intrinsic;

/*
intrinsic ReduceWeightedProjectiveSpace(IC::SeqEnum) -> SeqEnum
  {}
  dens:=[ Denominator(a) : a in IC ];
  p:=2;
  pd:=[ Valuation(a,p) : a in dens ];
  [ vs[1] - (vs[1] mod 2), vs[2] mod 4, vs[3] mod 6, vs[4] mod 10 ];
  end intrinsic
*/
