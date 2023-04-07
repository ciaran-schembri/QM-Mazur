
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


intrinsic UnitGroup(OmodN::AlgQuatOrdRes) -> GrpPerm, Map
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



intrinsic UnitGroup(O::AlgQuatOrd,N::RngIntElt) -> GrpPerm, Map
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
  UnitElements:=Set([ phi(x) : x in Setseq(Set(ONx)) ]);
  ZmodN:=ResidueClassRing(N);
 
  auts:=Set([ AutmuO(a) : a in Setseq(Set(Domain(AutmuO))) ]);

  enhancedimage:=Setseq(Set(CartesianProduct(auts, UnitElements)));

  semidirGL4xGL4seq:= [ <NormalizingElementToGL4modN(s[1],OmodN : basis:=basis),
     UnitGroupModNToGL4(s[2] : basis:=basis)> : s in enhancedimage ];
  semidirGL4seq:= [ a[1]*a[2] : a in semidirGL4xGL4seq ];
  semidirGL4:= sub< GL(4,ZmodN) |  semidirGL4seq >;
  assert #Set(semidirGL4) eq #semidirGL4xGL4seq;

  mapfromenhancedimage := map<  enhancedimage -> semidirGL4xGL4seq  |  
    s :-> <NormalizingElementToGL4modN(s[1],OmodN : basis:=basis), UnitGroupModNToGL4(s[2] : basis:=basis)>, 
    h :-> enhancedimage[Index(semidirGL4xGL4seq,h)]  >;
  maptogroup:= map< semidirGL4xGL4seq -> semidirGL4 |  a :-> a[1]*a[2], g :-> semidirGL4xGL4seq[Index(semidirGL4seq,g)]  >;

  EnhancedImageToGL4 := mapfromenhancedimage*maptogroup;

  return semidirGL4,EnhancedImageToGL4;
end intrinsic;


intrinsic EnhancedImageGL4(AutmuO::Map, O::AlgQuatOrd, N::RngIntElt) -> GrpMat
  {}

  OmodN:=quo(O,N);
  return EnhancedImageGL4(AutmuO, OmodN);
end intrinsic;


intrinsic FixedSubspace(H::GrpMat) -> Any 
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


/*
B<i,j,ij> := QuaternionAlgebra(Rationals(),-3,2);
BxmodQx:=QuaternionAlgebraModuloScalars(B);
O:=MaximalOrder(B);
mu := i*(2+j);
C2:=AbelianGroup([2]);
AutmuO:=map< C2 -> BxmodQx | [<C2!1,mu>,<C2!0,1>] >;
N:=4;
OmodN:=quo(O,N);

G,embed:=EnhancedImageGL4(AutmuO,O,N);
subs:=Subgroups(G);
for H in subs do
  Hgp:=H`subgroup;
  gens:=Generators(Hgp);
  gens_enhanced := [ Inverse(embed)(g) : g in gens ];
  order:=H`order;
  index:=Order(G)/order;
  fixedspace:=FixedSubspace(Hgp);

  printf "%o | %o | %o | %o \n", PrimaryAbelianInvariants(fixedspace), order,index,gens_enhanced;
end for;


*/





