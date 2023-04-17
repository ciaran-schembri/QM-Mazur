

B:=QuaternionAlgebra(6);
O:=MaximalOrder(B);
N:=3;

SetProfile(true);
for del in Divisors(Discriminant(B)) do
  if HasPolarizedElementOfDegree(O,del) then 
    tr,mu:=HasPolarizedElementOfDegree(O,del);
    rhocirc:=AllEnhancedSubgroups(O,mu,N: minimal:=true, verbose:=true, PQMtorsion:=true);
    G:=ProfileGraph();
    ProfilePrintByTotalTime(G : Max:=10);
  end if;
 end for;

del:=1;
tr,mu:=HasPolarizedElementOfDegree(O,del);
AutFull:=Aut(O,mu);
G,Gelts:=EnhancedImageGL4(AutFull,O,N);
H := rhocirc[8];
Hgp:=H`subgroup;
T := CosetTable(G,Hgp);
piH := CosetTableToRepresentation(G,T);
D:=Domain(AutFull);


B:=QuaternionAlgebra(6);
BxmodQx:=QuaternionAlgebraModuloScalars(B);
O:=MaximalOrder(B);
N:=3;

del:=1;
tr,mu:=HasPolarizedElementOfDegree(O,del);
AutFull:=Aut(O,mu);
OmodN:=quo(O,N);
ZmodN:=ResidueClassRing(N);


Fuchs:=FuchsianGroup(O);
FP,m1,m2:=Group(Fuchs);
FPgens:=Generators(FP);
FPfin:=[ <eval(Prune(Prune(Sprint(LHS(rel))))),eval(Sprint(LHS(rel))[#Sprint(LHS(rel))]) > : rel in Relations(FP) ];
Oxfin:=[ <m2(w[1]),w[2]> : w in FPfin ];
assert Sort([ w[2] : w in FPfin ]) eq Signature(Fuchs)[2];
[ Order(BxmodQx!(B!(a[1]*chi))) : a in Oxfin ];

AutmuOset:=[ AutFull(w) : w in Domain(AutFull)];
[ Norm(a`element) : a in AutmuOset ];
elliptic_elements_init:= [ <BxmodQx!w,BxmodQx!(B!x[1])> : w in AutmuOset, x in Oxfin ];
elliptic_elements:= [ a[1]*a[2] : a in elliptic_elements_init ];

normalizer_elliptic:=[ a : a in elliptic_elements | Type(Order(a)) eq RngIntElt ];
[Order(a) : a in normalizer_elliptic ];

Rx<x>:=PolynomialRing(Rationals());
E<u>:=NumberField(x^2+x-1);
v:=(1-u)/2;
Rm:=Order([1,u]);

chi,map:=Embed(Rm,O);
Norm(chi);
  cyc,Czeta,zeta:=IsCyclotomic(Em);
 

elliptic_inGL4:=[ ]

AutGO:=AutomorphismGroup(GO);
AutOmodqx := AutomorphismGroup(Omodqx);
theta := hom<AutmuO -> AutOmodqx | [AutOmodqx | thetaj,thetamuj]>;
AutmuinGL4:=sub< GL(4,ZmodN) | Setseq(Set([ a`GL4 : a in Gelts | a`enhanced[2] eq OmodN!(O!1) ])) >;
gensAutmuGL4:=Setseq(Generators(AutmuinGL4));


sigma:= [ piH(gamma) : gamma in gensAutmuGL4 ];
d := H`index;
// Riemann-Hurwitz formula
g2 := -2*d + &+[&+[ e[2]*(e[1]-1) : e in CycleStructure(sigma[i])] : i in [1..#sigma]];
g := (g2+2)/2;
g;



omegaO := Eltseq(K!S.2)[1] + Eltseq(K!S.2)[2]*(-1+i)/2;
O := QuaternionOrder([1,omegaO,j,omegaO*j]);
delta2 := mu;
delta4 := 1 - i + 1/2*j - 1/2*ij;
delta6 := Conjugate(delta2*delta4);
eps := 1+j;
assert IsScalar(delta2*delta4*delta6) and IsScalar(delta2^2) and IsScalar(delta4^4) and IsScalar(delta6^6);

g2 := iotaAutmuO(AutmuO![1,1]);
g4 := iotaOmodqx(matmodq(delta4/j))*iotaAutmuO(AutmuO![1,0]);
g6 := iotaOmodqx(matmodq(delta6/(mu*j)))*iotaAutmuO(AutmuO![0,1]);
eps := iotaOmodqx(matmodq(eps));
G1 := sub<G | [g2,g4,g6,eps]>;

Gamma1, m1 := Group(FuchsianGroup(O));


for Hrec in Subgroups(G) do
  // still have to think about a conjugate under Aut_mu(O)
  H := Hrec`subgroup;
  // if sub<Uq | [Determinant(GL4red(h))@@mq : h in Generators(H)]> ne Uq then continue; end if;    
  T := CosetTable(G,H);
  piH := CosetTableToRepresentation(G,T);
  sigma := [piH(g2),piH(g4),piH(g6)];
  fixsub := &meet[Kernel(GL4red(h)-1) : h in {Id(H)} join Generators(H)];
  abinvs := [];
  for i := 1 to Rank(fixsub) do
    j := 1; while (fixsub.i)[j] eq 0 do j +:= 1; end while;
    Append(~abinvs, q/Integers()!(fixsub.i)[j]);
  end for;
  // if 4 in abinvs then 
    print #H, Index(G,H), #(iotaOmodqx(Omodqx) meet H), #proj(H), Genus(sigma), abinvs;
  // end if;
end for;



intrinsic Genus(sigma::SeqEnum);
  {Compute genus from permutation triple}
  d := Degree(Parent(sigma[1]));
  // Riemann-Hurwitz formula
  g2 := -2*d + &+[&+[ e[2]*(e[1]-1) : e in CycleStructure(sigma[i])] : i in [1..3]];
  assert g2 mod 2 eq 0;
  g := (g2+2) div 2;
  return g;
end intrinsic;

