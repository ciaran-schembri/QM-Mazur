


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


intrinsic AtkinLehnerToAutmuO(a::AlgQuatElt,mu::AlgQuatElt,O::AlgQuatOrd) -> AlgQuatEnhElt 
  {}
  Ocirc:=EnhancedSemidirectProduct(O);
  AutFull,autmuOseq:=Aut(O,mu);

  assert a^2/Norm(a) in O;
  assert Norm(a) gt 0;
  for w in autmuOseq do 
    if IsSquare(Rationals()!Abs(Norm((w`element)^-1*a))) then
      w1:=w;
      tr,c:=IsSquare(Rationals()!Abs(Norm((w`element)^-1*a)));
      x:=(1/c)*((w`element)^-1)*a;
      assert x in O;
      assert Norm(x) in {1,-1};
      ell:=Ocirc!<w,O!x>;
      return ell;
    end if;
  end for;
end intrinsic;

intrinsic AtkinLehnerToAutmuO(a::AlgQuatElt,mu::AlgQuatOrdElt,O::AlgQuatOrd) -> AlgQuatEnhElt 
  {}
  return AtkinLehnerToAutmuO(a,QuaternionAlgebra(O)!mu,O);
end intrinsic;


intrinsic AtkinLehnerToAutmuO(a::AlgQuatOrdElt,mu::AlgQuatOrdElt,O::AlgQuatOrd) -> AlgQuatEnhElt 
  {}
  return AtkinLehnerToAutmuO(QuaternionAlgebra(O)!a,QuaternionAlgebra(O)!mu,O);
end intrinsic;

intrinsic AtkinLehnerToAutmuO(a::AlgQuatOrdElt,mu::AlgQuatElt,O::AlgQuatOrd) -> AlgQuatEnhElt 
  {}
  return AtkinLehnerToAutmuO(QuaternionAlgebra(O)!a,mu,O);
end intrinsic;



intrinsic NormalizerPlusGenerators(O::AlgQuatOrd) -> SeqEnum 
  {return generators of the positive norm elements which normalize O}
  if Discriminant(O) eq 6 then 
    B6<i6,j6>:=QuaternionAlgebra<Rationals() | -1,3 >;
    B:=QuaternionAlgebra(O);
    tr,map:=IsIsomorphic(B6,B : Isomorphism:=true);
    assert tr;
    B6elliptic_elts:=[ 3+3*i6+j6+i6*j6, 1+i6, 3*i6 + i6*j6];
    Oelliptic_elts:=[ O!map(a) : a in B6elliptic_elts ];
    assert Set([ Norm(a) : a in Oelliptic_elts ]) eq {2,6,12};
    return Oelliptic_elts;
  else 
    return "oops, not written for this discriminant yet";
  end if;
end intrinsic;


intrinsic NormalizerPlusGeneratorsEnhanced(O::AlgQuatOrd,mu::AlgQuatElt) -> Tup 
  {return generators of the positive norm elements which normalize O in the enhanced semidirect product}
  return < AtkinLehnerToAutmuO(a,mu,O) : a in NormalizerPlusGenerators(O) >;
end intrinsic;

intrinsic NormalizerPlusGeneratorsEnhanced(O::AlgQuatOrd,mu::AlgQuatOrdElt) -> Tup 
  {return generators of the positive norm elements which normalize O in the enhanced semidirect product}
  return < AtkinLehnerToAutmuO(a,mu,O) : a in NormalizerPlusGenerators(O) >;
end intrinsic;





intrinsic EnhancedEllipticElements(O::AlgQuatOrd,mu::AlgQuatElt) -> SeqEnum 
  {return the elliptic elements}
  if Discriminant(O) eq 6 then 
  return < AtkinLehnerToAutmuO(a,mu,O) : a in NormalizerPlusGenerators(O) >;
  else 
    return "oops, not written for this discriminant yet";
  end if;
end intrinsic;

intrinsic EnhancedEllipticElements(O::AlgQuatOrd,mu::AlgQuatOrdElt) -> SeqEnum
  {return the elliptic elements of the enhanced semidirect product}

  B:=QuaternionAlgebra(O);
  return EnhancedEllipticElements(O,B!mu);
end intrinsic;




