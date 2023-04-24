


intrinsic EnhancedGenus(sigma::SeqEnum) -> RngIntElt
  {Compute genus from permutation triple
   f:X -> Y. 2gX-2 = deg(f)*(2gY-2) + sum_x\inX (ex -1). 
   ex is the ramification degree of x. An element of S_n acts on sheets of the cover. 
  x is ramified if x is sent to another point under the action of an isotropy subgroup,
  i.e. the cycle type corresponding to x has length >1. The length is the ramification degree.}
  d := Degree(Parent(sigma[1]));
  // Riemann-Hurwitz formula
  rhs := -2*d + &+[ &+[ e[2]*(e[1]-1) : e in CycleStructure(sig) ] : sig in sigma ];
  assert rhs mod 2 eq 0;
  g := Integers()!((rhs+2)/2);
  return g;
end intrinsic;

intrinsic EnhancedCosetRepresentation(G::GrpMat,H::GrpMat) -> HomGrp
  {}
  if -H!1 notin H then 
    Hnew:=sub<G | H, -H!1 >;
  else 
    Hnew:=H;
  end if;
  T := CosetTable(G,Hnew);
  piH := CosetTableToRepresentation(G,T);
  return piH;
end intrinsic;


/*
intrinsic EnhancedEllipticElements(O::AlgQuatOrd,mu::AlgQuatElt) -> SeqEnum
  {return the elliptic elements of the enhanced semidirect product}
  Ocirc:=EnhancedSemidirectProduct(O);
  Fuchs:=FuchsianGroup(O);
  FP,m1,m2:=Group(Fuchs);
  FPgens:=Generators(FP);
  FPfin:=[ <eval(Prune(Prune(Sprint(LHS(rel))))),eval(Sprint(LHS(rel))[#Sprint(LHS(rel))]) > : rel in Relations(FP) ];
  Oxfin:=[ <m2(w[1]),w[2]> : w in FPfin ];
  assert Sort([ w[2] : w in FPfin ]) eq Signature(Fuchs)[2];
  //[ Order(BxmodQx!(B!(a[1]*chi))) : a in Oxfin ];

  AutFull,autmuOseq:=Aut(O,mu);

  elliptic_elements_init:=[];
  for x in Oxfin do
    for w in autmuOseq do 
      z:=x[1]+1; 
      if IsSquare(Rationals()!Abs(Norm((w`element)^-1*z))) then
        w1:=w;
        tr,c:=IsSquare(Rationals()!Abs(Norm((w`element)^-1*z)));
        assert (1/c)*((w`element)^-1)*z in O; 
        ell:=Ocirc!<w,O!((1/c)*((w`element)^-1)*z)>;
        Append(~elliptic_elements_init,ell);
      end if;
    end for;
  end for;

  elliptic_elements2:=[];

  for ell1,ell2 in elliptic_elements_init do 
    if ell1 ne ell2 and Type(Order(ell1*ell2)) eq RngIntElt then 
      Append(~elliptic_elements,ell1*ell2);
      //Append(~orders,Order(ell)); 
    end if;
  end for;

  for a in elliptic_elements2 do
    if exists(e){ b : b in elliptic_elements_init | a in [b^i : i in [1..Order(b)] ] } then 
      ;
    else 
      //a; Order(a);
      Append(~elliptic_elements_init,a);
    end if;
  end for;

  elliptic_elements:=[];
  orders:=[];
  for ell in elliptic_elements_init do 
    if Order(ell) notin orders then
      Append(~elliptic_elements,ell);
      Append(~orders,Order(ell));
    end if;
  end for;

  return elliptic_elements;

end intrinsic;
*/


intrinsic EnhancedEllipticElements(O::AlgQuatOrd,mu::AlgQuatElt) -> SeqEnum 
  {return the elliptic elements}
  if Discriminant(O) eq 6 then 
    B6<i6,j6>:=QuaternionAlgebra<Rationals() | -1,3 >;
    B:=QuaternionAlgebra(O);
    tr,map:=IsIsomorphic(B6,B : Isomorphism:=true);
    assert tr;
    B6elliptic_elts:=[ 3+3*i6+j6+i6*j6, 1+i6, 3*i6 + i6*j6];
    Oelliptic_elts:=[ O!map(a) : a in B6elliptic_elts ];
    assert Set([ Norm(a) : a in Oelliptic_elts ]) eq {2,6,12};
  else 
    return "oops, not written for this discriminant yet";
  end if;
end intrinsic;

intrinsic EnhancedEllipticElements(O::AlgQuatOrd,mu::AlgQuatOrdElt) -> SeqEnum
  {return the elliptic elements of the enhanced semidirect product}

  B:=QuaternionAlgebra(O);
  return EnhancedEllipticElements(O,B!mu);
end intrinsic;

