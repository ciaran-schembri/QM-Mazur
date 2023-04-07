B<i,j,ij> := QuaternionAlgebra(Rationals(),-3,2);
BxmodQx:=QuaternionAlgebraModuloScalars(B);
O:=MaximalOrder(B);
mu := i*(2+j);
C2:=AbelianGroup([2]);
AutmuO:=map< C2 -> BxmodQx | [<C2!1,mu>,<C2!0,1>] >;
N:=6;
OmodN:=quo(O,N);

G,embed:=EnhancedImageGL4(AutmuO,O,N);
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

