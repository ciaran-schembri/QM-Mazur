
declare type AlgQuatOrdRes[AlgQuatOrdResElt];

declare attributes AlgQuatOrdRes :
  quaternionorder,
  quaternionideal;

declare attributes AlgQuatOrdResElt:
  element,
  parent;

declare type AlgQuatProj[AlgQuatProjElt];

declare attributes AlgQuatProj :
  quaternionalgebra;

declare attributes AlgQuatProjElt :
  element,
  parent;

intrinsic OmodNElement(OmodN::AlgQuatOrdRes, x::AlgQuatOrdElt) -> AlgQuatOrdResElt
  {Construct an element of the OmodN whose underlying element is x in O}
  elt := New(AlgQuatOrdResElt);
  elt`element := x;
  elt`parent := OmodN;
  
  return OmodN!elt;
end intrinsic;

intrinsic ElementModuloScalars(BxmodFx::AlgQuatProj, x::AlgQuatElt) -> AlgQuatProjElt
  {Construct an element of B^x/F^x whose underlying element is x in B}
  elt := New(AlgQuatProjElt);
  elt`element := x;
  elt`parent := BxmodFx;
  
  return elt;
end intrinsic;

intrinsic 'eq'(x::AlgQuatOrdResElt,y::AlgQuatOrdResElt) -> BoolElt 
  {Decide if x equals y in OmodN}
  assert Parent(x) eq Parent(y);
  OmodN:=Parent(x);
  N:=OmodN`quaternionideal;
  O:=OmodN`quaternionorder;
  xO:=x`element;
  yO:=y`element;
  return xO-yO in N*O;
end intrinsic;

intrinsic 'eq'(x::AlgQuatProjElt,y::AlgQuatProjElt) -> BoolElt 
  {Decide if x equals y in OmodN}
  assert Parent(x) eq Parent(y);
  //BxmodFx:=Parent(x);
  x0:=x`element;
  y0:=y`element;
  assert x0*y0 ne 0;
  return IsScalar(x0/y0);
end intrinsic;

intrinsic 'eq'(OmodN1::AlgQuatOrdRes,OmodN2::AlgQuatOrdRes) -> BoolElt 
  {Decide if OmodN1 equals OmodN2}
 
  O1:=OmodN1`quaternionorder;
  O2:=OmodN2`quaternionorder;

  N1:=OmodN1`quaternionideal;
  N2:=OmodN2`quaternionideal;

  if O1 eq O2 and N1 eq N2 then 
    return true;
  else 
    return false;
  end if;
end intrinsic;


intrinsic 'eq'(BxmodFx1::AlgQuatProj,BxmodFx2::AlgQuatProj) -> BoolElt 
  {Decide if BxmodFx1 equals BxmodFx2}
 
  B1:=BxmodFx1`quaternionalgebra;
  B2:=BxmodFx2`quaternionalgebra;

  if B1 eq B2 then 
    return true;
  else 
    return false;
  end if;
end intrinsic;

intrinsic '*'(x::AlgQuatOrdResElt,y::AlgQuatOrdResElt) -> AlgQuatOrdResElt 
  {compute x*y in OmodN}
  assert Parent(x) eq Parent(y);
  OmodN:=Parent(x);
  N:=OmodN`quaternionideal;
  O:=OmodN`quaternionorder;
  xO:=x`element;
  yO:=y`element;

  ON,piN:=quo(O,N);
  return OmodN!(piN(xO*yO));
end intrinsic;

intrinsic '*'(x::AlgQuatProjElt,y::AlgQuatProjElt) -> AlgQuatProjElt 
  {compute x*y in B^x/F^x}
  assert Parent(x) eq Parent(y);
  BxmodFx:=Parent(x);

  xO:=x`element;
  yO:=y`element;

  return BxmodFx!ElementModuloScalars(BxmodFx,xO*yO);
end intrinsic;

intrinsic '^'(x::AlgQuatProjElt,y::RngIntElt) -> AlgQuatProjElt 
  {compute x*y in B^x/F^x}
  BxmodFx:=Parent(x);

  xO:=x`element;

  return BxmodFx!ElementModuloScalars(BxmodFx,xO^y);
end intrinsic;
  
intrinsic Parent(elt::AlgQuatOrdResElt) -> AlqQuatOrdRes
  {.}
  return elt`parent;    
end intrinsic;

intrinsic Parent(elt::AlgQuatProjElt) -> AlqQuatProj
  {.}
  return elt`parent;    
end intrinsic;

intrinsic quo(O::AlgQuatOrd, N::RngIntElt) -> AlgQuatOrdRes
  {.}
  M := New(AlgQuatOrdRes);
  M`quaternionideal := N;
  M`quaternionorder := O;
  projection := map < O -> M | x :-> OmodNElement(M,x) >;
  return M, projection;
end intrinsic;

intrinsic QuaternionAlgebraModuloScalars(B::AlgQuat) -> AlgQuatProj 
  {Create B^x/F^x}
  BxmodFx:=New(AlgQuatProj);
  BxmodFx`quaternionalgebra := B;
  return BxmodFx;
end intrinsic;

intrinsic IsCoercible(OmodN::AlgQuatOrdRes, x::Any) -> BoolElt, .
{.}
  if Type(x) eq AlgQuatOrdResElt then
    if Parent(x) eq OmodN then
      return true, x;
    else
      return false, "Illegal Coercion";
    end if;
  elif Type(x) eq AlgQuatOrdElt then
    if Parent(x) eq OmodN`quaternionorder then 
      return true,OmodNElement(OmodN,x);
    else 
      return false, "Illegal Coercion";
    end if;
  else
    return false, "Illegal Coercion";
  end if;
end intrinsic;


intrinsic IsCoercible(BxmodFx::AlgQuatProj, x::Any) -> BoolElt, .
{.}
 if Type(x) eq AlgQuatProjElt then
    if Parent(x) eq BxmodFx then
      return true, x;
    else
      return false, "Illegal Coercion";
    end if;
  elif Type(x) eq AlgQuatElt then 
    if Parent(x) eq BxmodFx`quaternionalgebra then 
      return true, ElementModuloScalars(BxmodFx,x);
    else   
      return false, "Illegal Coercion";   
    end if;
  else
    return false, "Illegal Coercion";
  end if;

end intrinsic;



intrinsic IsUnit(x::AlgQuatOrdResElt) -> BoolElt
  {return whether x \in O/N is a unit}
  x0:=x`element;
  OmodN:=Parent(x);
  N:=OmodN`quaternionideal;
  nm:=Norm(x0);
  ZmodN:=ResidueClassRing(N);

  if IsUnit(ZmodN!nm) then 
    return true;
  else 
    return false;
  end if;
end intrinsic;

intrinsic Set(OmodN::AlgQuatOrdRes) -> Set 
  {return the set of elements O/N}

  O:=OmodN`quaternionorder;
  N:=OmodN`quaternionideal;
  ON,piN:=quo(O,N);
  basis:=Basis(O);
  set:={  O!(a*basis[1] + b*basis[2] + c*basis[3] + d*basis[4]) : a,b,c,d in [0..N-1]  };
  return { OmodN!piN(x) : x in set };
end intrinsic;


/*function RightRegularRepresentation(U)
  G := Sym(#U);
  Useq:=Setseq(U);
  permrep:=map< U -> G | u :-> G![ Index(Useq,x*u) : x in Useq ] >;
  perms := [ permrep(u) : u in Useq ];
  S:=sub< G | perms >;
  //assert #S eq #U;
  return S, permrep;
end function;*/
/*  [ G | [ Index(Useq,x*u) : x in Useq ] : u in Useq ];
    H := sub< G | perms[gens] >;
    if not IsEmpty(K) then
        K := sub< H | perms[K] >;
        m := RegularRepresentation(H,K);
        perms:= perms @ m;
        H := Codomain(m);
    end if;
    assert #H eq #U;
    return H, map< H -> U | h:-> U[Index(perms, h)], u:-> perms[Index(U,u)] >;
end function;*/


intrinsic UnitGroup(OmodN::AlgQuatOrdRes) -> GrpMat, Map
  {return (O/N)^x as a permutation group G, the second value is the isomorphism G ->(O/N)^x}
  //Need to make this much more efficient.

  O:=OmodN`quaternionorder;
  N:=OmodN`quaternionideal;
  units := { x : x in Set(OmodN) | IsUnit(x) };
  Useq:=Setseq(units);

  unitsinGL4:= [ UnitGroupModNToGL4(x) : x in Useq ];

  ZmodN:=ResidueClassRing(N);

  subONx:=sub< GL(4,ZmodN) | unitsinGL4 >;
  assert #Set(subONx) eq #units;

  phi:=map< subONx -> Useq | s :-> Useq[Index(unitsinGL4,s)], x :-> UnitGroupModNToGL4(x) >;

  /*G := Sym(#units);
  permrep:=map< units -> G | u :-> G![ Index(Useq,x*u) : x in Useq ] >;
  perms := [ permrep(u) : u in Useq ];

  k:=1;
  S:=sub< G | permrep(Useq[k]) >;
  while #S lt #units do 
    k:=k+1;
    S:=sub<G | permrep([ Useq[l] : l in [1..k]]) >;
  end while;
  assert #S eq #units;

  //assert group is UU is correct.
  phi := map < S -> units | u :-> Useq[Index(perms,u)], x :-> permrep(x) >;
  //phi_inv := map< units -> S | x :-> permrep(x) >;*/
  return subONx,phi;
end intrinsic;



intrinsic UnitGroup(O::AlgQuatOrd,N::RngIntElt) -> GrpMat, Map
  {return (O/N)^x as a permutation group G, the second value is the isomorphism G ->(O/N)^x}

  return UnitGroup(quo(O,N));
end intrinsic;



intrinsic ElementToAutomorphismModN(a::AlgQuatElt, OmodN::AlgQuatOrdRes) -> GrpAutoElt
  {a in B^x becomes an automorphism of (O/N)^x by considering the map a |-> (x|-> a^-1xa) 
  as long as a \in N_B^x(O). We apply this to (O/N)^x as a permutation group.}

  O:=OmodN`quaternionorder;
  N:=OmodN`quaternionideal;
  ON,piN:=quo(O,N);
  ONx,phi,phi_inv:=UnitGroup(OmodN);
  aut:=hom< ONx -> ONx | x :-> phi_inv(piN(a^(-1)*(phi(x)`element)*a)) >;
  return AutomorphismGroup(ONx)!aut;
end intrinsic;

intrinsic ElementToAutomorphismModN(a::AlgQuatElt, O::AlgQuatOrd, N::RngIntElt) -> GrpAutoElt
  {a in B^x becomes an automorphism of (O/N)^x by considering the map a |-> (x|-> a^-1xa) 
  as long as a \in N_B^x(O). We apply this to (O/N)^x as a permutation group.}

  return ElementToAutomorphismModN(a,quo(O,N));
end intrinsic;


intrinsic AutomorphismsModN(S::{ AlgQuatProjElt }, OmodN::AlgQuatOrdRes) -> Map
  {Given a subset of Aut(O) input as a finite subset S of B^x, create the map
   theta : S -> Aut((O/N)^x)}
  
  ONx,phi:=UnitGroup(OmodN);
  return map< S -> AutomorphismGroup(ONx) | s:->ElementToAutomorphismModN(s,OmodN) >;
end intrinsic;

intrinsic AutomorphismsModN(S::{ AlgQuatProjElt }, O::AlgQuatOrd, N::RngIntElt) -> Map
  {Given a subset of Aut(O) input as a finite subset S of B^x, create the map
  theta : S -> Aut((O/N)^x)}
  
  OmodN:=quo(O,N);
  return AutomorphismsModN(S,OmodN);
end intrinsic;


intrinsic MapIsHomomorphism(AutmuO::. : injective:=true) -> BoolElt
  {Check whether the map AutmuO : C -> B^x/Q^x is an injective homomorphism}
  for x,y in Domain(AutmuO) do 
    if not(AutmuO(x*y) eq AutmuO(x)*AutmuO(y)) then 
      return false;
    end if;
    if injective eq true then 
      if ((AutmuO(x) eq AutmuO(y)) and (x ne y)) then 
        return false;
      end if;
    end if;
  end for;
  return true;
end intrinsic 




intrinsic NormalizingElementToGL4modN(x::AlgQuatElt,OmodN::AlgQuatOrdRes : basis:=[]) -> GrpMatElt 
  {O is an order over R. For an element g \in N_Bx(O) the map phi_g : b |--> g^-1bg
  is R-linear hence [g] is an element of M_4(R) after fixing a basis
  this function computes [g] and also returns the R-basis of O.}

    O:=OmodN`quaternionorder;
    N:=OmodN`quaternionideal;
    if basis eq [] then 
      basis:=Basis(O);
    end if;
    R:=BaseRing(O);
    assert R eq Integers();
    M4R:=MatrixAlgebra(R,4);
    ZmodN:=ResidueClassRing(N);

    x_map:=Transpose(M4R![ Eltseq(O!(x^(-1)*b*x)) : b in basis ]);
    assert Determinant(x_map) eq 1;

    return GL(4,ZmodN)!x_map, basis;
end intrinsic;


intrinsic NormalizingElementToGL4modN(x::AlgQuatProjElt,OmodN::AlgQuatOrdRes : basis:=[]) -> GrpMatElt 
  {O is an order over R. For an element g \in N_Bx(O) the map phi_g : b |--> g^-1bg
  is R-linear hence [g] is an element of M_4(R) after fixing a basis
  this function computes [g] and also returns the R-basis of O.}

  return NormalizingElementToGL4modN(x`element, OmodN : basis:=basis );
end intrinsic;


intrinsic NormalizingElementToGL4modN(x::AlgQuatElt,O::AlgQuatOrd, N::RngIntElt : basis:=[]) -> GrpMatElt 
  {O is an order over R. For an element g \in N_Bx(O) the map phi_g : b |--> g^-1bg
  is R-linear hence [g] is an element of M_4(R) after fixing a basis
  this function computes [g] and also returns the R-basis of O.}

  OmodN:=quo(O,N);
  return NormalizingElementToGL4modN(x,OmodN : basis:=basis);
end intrinsic;
  

intrinsic NormalizingElementToGL4modN(x::AlgQuatProjElt,O::AlgQuatOrd, N::RngIntElt : basis:=[]) -> GrpMatElt 
  {O is an order over R. For an element g \in N_Bx(O) the map phi_g : b |--> g^-1bg
  is R-linear hence [g] is an element of M_4(R) after fixing a basis
  this function computes [g] and also returns the R-basis of O.}

  OmodN:=quo(O,N);
  return NormalizingElementToGL4modN(x,OmodN : basis:=basis);
end intrinsic;


intrinsic UnitGroupModNToGL4(x::AlgQuatOrdResElt : basis:=[]) -> GrpMatElt 
  {O is an order over R, this returns a matrix [lambda_g] wrt to a basis
  which is the right regular representation
  lambda_g : g --> b*g where g \in GL_1(O)}

  OmodN:=Parent(x);
  x0:=x`element;
  O:=OmodN`quaternionorder;
  N:=OmodN`quaternionideal;
  if basis eq [] then 
    basis:=Basis(O);
  end if;
  R:=BaseRing(O);
  assert R eq Integers();
  M4R:=MatrixAlgebra(R,4);
  ZmodN:=ResidueClassRing(N);


  x_map:=Transpose(M4R![ Eltseq(O!(b*x0)) : b in basis ]);
  assert ZmodN!Determinant(x_map) ne 0;
  return GL(4,ZmodN)!x_map;
end intrinsic;
 


intrinsic EnhancedImagePermutation(AutmuO::Map,OmodN::AlgQuatOrdRes) -> Grp 
  {AutmuO is a map from a finite group C -> B^x, which is isomorphic onto the image in B^x/Q^x. 
  We create the semidirect product of ONx by AutmuO, using AutomorphismModN as the map
  theta: AutmuO -> Aut(ONx)}

  assert MapIsHomomorphism(AutmuO : injective:=true);
  H:=Domain(AutmuO);
  ONx,phi := UnitGroup(OmodN);
  AutONx:=AutomorphismGroup(ONx);
  theta:=hom< H -> AutONx | h :-> ElementToAutomorphismModN(AutmuO(h),OmodN) >;

  rho_circ,m1,m2,m3:=SemidirectProduct(ONx,H,theta);
  return rho_circ,m1,m2,m3;
end intrinsic;

intrinsic EnhancedImagePermutation(AutmuO::.,O::AlgQuatOrd, N::RngIntElt) -> Grp 
  {AutmuO is a map from a finite group C -> B^x, which is isomorphic onto the image in B^x/Q^x. 
  We create the semidirect product of ONx by AutmuO, using AutomorphismModN as the map
  theta: AutmuO -> Aut(ONx)}

  OmodN:=quo(O,N);
  return EnhancedImagePermutation(AutmuO,OmodN);
end intrinsic;


intrinsic EnhancedImageGL4(AutmuO::Map, OmodN::AlgQuatOrdRes) -> GrpMat
  {}

  O:=OmodN`quaternionorder;
  N:=OmodN`quaternionideal;
  basis:=Basis(O);
  assert MapIsHomomorphism(AutmuO : injective:=true);
  H:=Domain(AutmuO);
  ONx,phi:= UnitGroup(OmodN);
  UnitElements:=[ phi(x) : x in ONx ];
  ZmodN:=ResidueClassRing(N);
 
  auts:=[ AutmuO(a) : a in Domain(AutmuO) ];

  enhancedimage_cartesian:=[ c: c in CartesianProduct(auts, UnitElements) ];

  RF := recformat< n : Integers(),
  enhanced,
  GL4xGL4,
  GL4
  >
  ;

  enhancedimage:=[];
  for elt in enhancedimage_cartesian do 
    s := rec< RF | >;
    s`enhanced:=elt;
    s`GL4xGL4:=<NormalizingElementToGL4modN(elt[1],OmodN : basis:=basis), UnitGroupModNToGL4(elt[2] : basis:=basis)>;
    s`GL4:=s`GL4xGL4[1]*s`GL4xGL4[2];
    Append(~enhancedimage,s);
  end for;

 /* semidirGL4xGL4seq:= [ <NormalizingElementToGL4modN(s[1],OmodN : basis:=basis),
     UnitGroupModNToGL4(s[2] : basis:=basis)> : s in enhancedimage ];
  semidirGL4seq:= [ a[1]*a[2] : a in semidirGL4xGL4seq ];
  semidirGL4:= sub< GL(4,ZmodN) |  semidirGL4seq >;
  assert #Set(semidirGL4) eq #semidirGL4xGL4seq;

  mapfromenhancedimage := map<  enhancedimage -> semidirGL4xGL4seq  |  
    s :-> <NormalizingElementToGL4modN(s[1],OmodN : basis:=basis), UnitGroupModNToGL4(s[2] : basis:=basis)>, 
    h :-> enhancedimage[Index(semidirGL4xGL4seq,h)]  >;
  maptogroup:= map< semidirGL4xGL4seq -> semidirGL4 |  a :-> a[1]*a[2], g :-> semidirGL4xGL4seq[Index(semidirGL4seq,g)]  >;
  EnhancedImageToGL4 := mapfromenhancedimage*maptogroup;
  */
  semidirGL4:= sub< GL(4,ZmodN) |  [ x`GL4 : x in enhancedimage ] >;

  return semidirGL4,enhancedimage;
end intrinsic;


intrinsic EnhancedImageGL4(AutmuO::Map, O::AlgQuatOrd, N::RngIntElt) -> GrpMat
  {}

  OmodN:=quo(O,N);
  return EnhancedImageGL4(AutmuO, OmodN);
end intrinsic;


intrinsic FixedSubspace(H::GrpMat) -> GrpAb 
  {}
  N:=#BaseRing(H);
  ZmodN:=ResidueClassRing(N);
  ZmodN4:= [ Matrix(ZmodN,1,4,[a,b,c,d]) : a,b,c,d in ZmodN ];
  fixed_vectors:= [ v : v in ZmodN4 | forall(u){ g : g in H | v*g eq v } ];
  A:=AbelianGroup([N,N,N,N]);
  elts:=[ Eltseq(v) : v in fixed_vectors ];
  fixedtors:=sub<A | elts >; 
  //fixedsub:=sub< AbelianGroup([N,N,N,N]) | 
  return fixedtors; //PrimaryAbelianInvariants(fixedsub);
end intrinsic;


intrinsic HasPolarizedElementOfDegree(O::AlgQuatOrd,d::RngIntElt) -> BoolElt, AlgQuatOrdElt 
  {return an element mu of O such that mu^2 + d*disc(O) = 0 if it exists.}
  disc:=Discriminant(O);
  Rx<x>:=PolynomialRing(Rationals());
  Em<v>:=NumberField(x^2+d*disc);
  if IsSplittingField(Em,QuaternionAlgebra(O)) then 
    cyc,Czeta,zeta:=IsCyclotomic(Em);
    if cyc eq true then      
      Zzeta:=Integers(Czeta);
      z:=Zzeta.2;
      zO,emb:=Embed(Zzeta,O);
      if CyclotomicOrder(Czeta) eq 4 then 
        assert zO^4 eq 1;
        tr,c:=IsSquare(Integers()!d*disc);
        mu:=c*zO;
        assert mu^2+d*disc eq 0;
        assert mu in O;
        return true, O!mu;
      else 
        assert zO^3 eq 1;
        tr,c:=IsSquare(Integers()!d*disc/3);
        mu:=(2*zO+1)*c;
        assert mu^2+d*disc eq 0;
        assert mu in O;
        return true,O!mu;
      end if;
    else 
      Rm:=Order([1,v]);
      mu,emb:=Embed(Rm,O);
      assert mu^2+d*disc eq 0;
      return true, O!mu;
    end if;
  else 
    return false;
  end if;
end intrinsic;


intrinsic DegreeOfPolarizedElement(O::AlgQuatOrd,mu:AlgQuatOrdElt) -> RngIntElt
  {degree of mu}
  tr,nmu:= IsScalar(mu^2);
  disc:=Discriminant(O);
  del:=-nmu/disc;
  assert IsCoercible(Integers(),del);
  assert IsSquarefree(Integers()!del);
  return Integers()!del;
end intrinsic;

intrinsic IsTwisting(O::AlgQuatOrd,mu::AlgQuatOrdElt) -> BoolElt
  {(O,mu) is twisting (of degree del = -mu^2/disc(O)) if there exists chi in O and N_Bx(O)
   such that chi^2 = m, m|Disc(O) and mu*chi = -chi*mu }

  tr,nmu:= IsScalar(mu^2);
  disc:=Discriminant(O);
  del:=DegreeOfPolarizedElement(O,mu);
  B:=QuaternionAlgebra(O);
  ram:=Divisors(disc) cat [ -1*m : m in Divisors(disc) ];

  for m in ram do
  Bram<i1,j1>:=QuaternionAlgebra< Rationals() | -disc*del, m>;
  tr,isom:=IsIsomorphic(Bram, B : Isomorphism:=true);
    if tr eq true then

      chi:=isom(j1);
      if mu*chi eq -chi*mu and chi in O then
        twisted_basis:=[1,mu,chi,mu*chi];
        Omuchi:=QuaternionOrder(Integers(), twisted_basis);   
        assert IsIsomorphic(O,Omuchi);
        return true, m,twisted_basis;
      end if;
    end if;
  end for;
  return false;
end intrinsic;



intrinsic Aut(O::AlgQuatOrd,mu::AlgQuatOrdElt) -> Any
  {}

  assert IsScalar(mu^2);
  tr,eta:=IsScalar(mu^2);
  disc:=Discriminant(O);
  Rx<x>:=PolynomialRing(Rationals());
  Em<v>:=NumberField(x^2-eta);
  //Rm:=Order([1,v]);
  cyc,Czeta,zeta:=IsCyclotomic(Em);
  //Zzeta:=Integers(Czeta);

  B:=QuaternionAlgebra(O);
  BxmodQx:=QuaternionAlgebraModuloScalars(B);

  if cyc eq true then  
    sqeta,c:=SquarefreeFactorization(Integers()!eta);
    assert sqeta in {-1,-3};
    if sqeta eq -1 then 
      cyc_order:=4;
      zeta_n := mu/c;
    elif sqeta eq -3 then 
      cyc_order:=6;
      zeta_n := ((mu/c)+1)/2;
    end if;
    a:=B!zeta_n+1;
  else 
    cyc_order:=2;
    a:=B!mu;
  end if;

  if IsTwisting(O,mu) then
    tr,m,twisted_basis:=IsTwisting(O,mu);
    b:=B!(twisted_basis[3]);
    if cyc eq true then 
      Dn:=DihedralGroup(cyc_order);
    else 
      Dn:=Group("C2^2");
    end if;
    Dngens:=Generators(Dn);
    assert #Dngens eq 2;
    assert Order(Dn.1) eq #Dn/2;
    assert Order(Dn.2) eq 2;
    elts:= [ <Dn.1^l*Dn.2^k, BxmodQx!(a^l*b^k)> : l in [0..cyc_order-1], k in [0..1] ];
    grp_map:=map< Dn -> BxmodQx | elts >;
  else 
    if cyc eq true then 
      Cn:=CyclicGroup(cyc_order);
    else 
      Cn:=CyclicGroup(2);
    end if;
    elts:= [ <Cn.1^k,BxmodQx!(a^k)> : k in [0..#Cn-1] ];
    grp_map:=map< Cn -> BxmodQx | elts >;
  end if;
 
  return grp_map;
end intrinsic;


intrinsic AllEnhancedSubgroups(O::AlgQuatOrd,mu::AlgQuatOrdElt,N::RngIntElt : minimal:=true,onlypotential:=true,verbose:=true) -> Any
  {}
  B:=QuaternionAlgebra(O);
  BxmodQx:=QuaternionAlgebraModuloScalars(B);
  OmodN:=quo(O,N);

  //mu:=PolarizedElementOfDegree(O,1);
  AutFull:=Aut(O,mu);
  assert MapIsHomomorphism(AutFull : injective:=true);

  RF := recformat< n : Integers(),
    subgroup,
    order,
    index,
    fixedsubspace,
    generators,
    split,
    endomorphism_representation,
    atkin_lehners
    >
    ;

  G,Gelts:=EnhancedImageGL4(AutFull,O,N);
  //components:=Components(embed);
  ZmodN:=ResidueClassRing(N);
  subs:=Subgroups(G);
  Autmuimage:=[AutFull(c) : c in Domain(AutFull) ];

  minimal_subs_init:=<>;

  for H in subs do
    Hgp:=H`subgroup;
    //Hgpset:= Set(Hgp);
    fixedspace:=FixedSubspace(Hgp);
    //if not(exists(e){ N : N in minimal_subs | Hgpset subset N and fixedspace eq }) then
      //Append(~minimal_subs,Hgpset);
    gens:=Generators(Hgp);
    gens_enhanced := [ g`enhanced : g in Gelts | g`GL4 in gens ];
    enhanced_set:= [ g`enhanced : g in Gelts | g`GL4 in Hgp ];

    for w in Autmuimage do 
      if <w,OmodN!(O!1)> in enhanced_set then 
        is_split := true;
      else 
        is_split := false;
      end if;
    end for;

    order:=H`order;
    index:=Order(G)/order;

    //endomorphism_image_set:=Set([ h[1] : h in enhanced_set ]);
    rho_end:=sub< GL(4,ZmodN) | Setseq(Set([ (g`GL4xGL4)[1] : g in Gelts | g`GL4 in Hgp ])) >;

    s := rec< RF | >;
    s`subgroup:=Hgp;
    s`order:=order;
    s`index:=index;
    s`fixedsubspace:=PrimaryAbelianInvariants(fixedspace);
    s`generators:=gens;
    s`split:=is_split;
    s`endomorphism_representation:=rho_end;
    s`atkin_lehners:=Sort([ SquarefreeFactorization(Integers()!Norm(x`element)) : x in Set([ (y`enhanced)[1] : y in Gelts | y`GL4 in Hgp  ]) ]);

    if onlypotential then 
      if #rho_end ne 1 then 
        Append(~minimal_subs_init,s);
      end if;
    else 
      Append(~minimal_subs_init,s);
    end if;
  end for;

  if minimal eq false then 
    return minimal_subs_init;
  else 
    minimal_subs:=<>;
    for s in minimal_subs_init do  
      F:=s`subgroup;
      tors:=s`fixedsubspace;
      endorep:=s`endomorphism_representation;
      AL:=s`atkin_lehners;
      if exists(e){ N : N in minimal_subs_init | F subset N`subgroup and 
        tors eq N`fixedsubspace and F ne N`subgroup 
         and AL eq N`atkin_lehners } then 
        ;
      else 
        Append(~minimal_subs,s);
      end if;
    end for;
    if verbose eq true then
      printf "Quaternion order of discriminant %o\n", Discriminant(O);
      printf "Polarized Element \\mu=%o of degree %o and norm %o\n", mu, DegreeOfPolarizedElement(O,mu),Norm(mu);
      for s in minimal_subs do 
        printf "%o | %o | %o | %o | %o | %o \n",  s`index, s`order, s`split, s`fixedsubspace, GroupName(s`endomorphism_representation), s`atkin_lehners;
      end for;
    end if;
    return minimal_subs;
  end if;

end intrinsic;


intrinsic AllEnhancedSubgroups(B::AlgQuat,mu::AlgQuatOrdElt,N::RngIntElt : minimal:=true,verbose:=true, onlypotential:=true) -> Any
  {}
  return AllEnhancedSubgroups(MaximalOrder(B),mu,N : minimal:=minimal, verbose:=verbose, onlypotential:=onlypotential);
end intrinsic;

intrinsic AllEnhancedSubgroups(O::AlgQuatOrd,del::RngIntElt,N::RngIntElt : minimal:=true,verbose:=true, onlypotential:=true) -> Any
  {}
  tr,mu:=HasPolarizedElementOfDegree(O,del);
  return AllEnhancedSubgroups(O,mu,N : minimal:=minimal, verbose:=verbose, onlypotential:=onlypotential);
end intrinsic;

intrinsic AllEnhancedSubgroups(B::AlgQuat,del::RngIntElt,N::RngIntElt : minimal:=true,verbose:=true, onlypotential:=true) -> Any
  {}
  O:=MaximalOrder(B);
  tr,mu:=HasPolarizedElementOfDegree(O,del);
  return AllEnhancedSubgroups(O,mu,N : minimal:=minimal, verbose:=verbose, onlypotential:=onlypotential);
end intrinsic;

intrinsic AllEnhancedSubgroups(D::RngIntElt,del::RngIntElt,N::RngIntElt : minimal:=true,verbose:=true, onlypotential:=true) -> Any
  {}
  B:=QuaternionAlgebra(D);
  O:=MaximalOrder(B);
  tr,mu:=HasPolarizedElementOfDegree(O,del);
  return AllEnhancedSubgroups(O,mu,N : minimal:=minimal, verbose:=verbose, onlypotential:=onlypotential);
end intrinsic;




intrinsic Print(elt::AlgQuatOrdResElt)
{.}
  printf "%o", elt`element;
end intrinsic;

intrinsic Print(OmodN::AlgQuatOrdRes)
{.}
  printf "Quotient of %o by %o", OmodN`quaternionorder, OmodN`quaternionideal;
end intrinsic;

intrinsic Print(elt::AlgQuatProjElt)
{.}
  printf "%o", elt`element;
end intrinsic;

intrinsic Print(BxmodFx::AlgQuatProj)
{.}
  printf "Quotient by scalars of %o", BxmodFx`quaternionalgebra;
end intrinsic;





