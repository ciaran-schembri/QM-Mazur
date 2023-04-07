B<i,j,ij> := QuaternionAlgebra(Rationals(),-3,2);
BxmodQx:=QuaternionAlgebraModuloScalars(B);
O:=MaximalOrder(B);
mu := i*(2+j);
C2:=CyclicGroup(2);
//AutmuO:=map< C2 -> BxmodQx | [<C2.1,mu>,<C2!1,1>] >;
N:=4;
OmodN:=quo(O,N);

mu:=PolarizedElementOfDegree(O,1);
AutFull:=Aut(O,mu);
assert MapIsHomomorphism(AutFull : injective:=true);

C:=Domain(AutFull);
subC:=Subgroups(C);
for S in Subgroups(C) do 
  s:=S`subgroup;
  s;
  [ SquarefreeFactorization(Integers()!Norm(AutFull(a)`element)) : a in s ];
end for;

AutmuO:=map< subC[2]`subgroup -> BxmodQx | a :-> AutFull(a) >;
assert MapIsHomomorphism(AutmuO : injective:=true);

SetProfile(true);

G,embed:=EnhancedImageGL4(AutFull,O,N);
T:=ProfileGraph();
ProfilePrintByTotalTime(T : Max:=10);

subs:=Subgroups(G);
for H in subs do
  Hgp:=H`subgroup;
  gens:=Generators(Hgp);
  gens_enhanced := [ Inverse(embed)(g) : g in gens ];
  enhanced_set:= [ Inverse(embed)(g) : g in Setseq(Set(Hgp)) ];
  is_split:=<mu,OmodN!(O!1)> in enhanced_set;
  order:=H`order;
  index:=Order(G)/order;
  fixedspace:=FixedSubspace(Hgp);

  printf "%o | %o | %o | %o \n", PrimaryAbelianInvariants(fixedspace), order,index, is_split;
end for;

