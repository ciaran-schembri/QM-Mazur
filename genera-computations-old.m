

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

////////////////////////////////////

B:=QuaternionAlgebra(6);
BxmodQx:=QuaternionAlgebraModuloScalars(B);
O:=MaximalOrder(B);
Ocirc:=EnhancedSemidirectProduct(O);
N:=3;

del:=1;
tr,mu:=HasPolarizedElementOfDegree(O,del);
rhocirc:=AllEnhancedSubgroups(O,mu,N: minimal:=true, verbose:=true, PQMtorsion:=true);

AutFull,autmuOseq:=Aut(O,mu);
[ Norm(a`element) : a in autmuOseq ];
OmodN:=quo(O,N);
ZmodN:=ResidueClassRing(N);


Qz6<z6>:=CyclotomicField(6);
Zz6:=Order([1,z6]);
c6:=Embed(Zz6,O);
Pc6:=1+c6;
assert IsScalar(Pc6^6);
chi:=autmuOseq[2];
chimu:=autmuOseq[4];
mub:=autmuOseq[3];
assert (chi`element)^2 eq 3;
a6:=O!(((chi`element)^-1)*Pc6);
//O!(((mub`element)^-1)*Pc6);
//O!(((chimu`element)^-1)*Pc6);
assert Norm(a6) eq -1;
//d6:=Ocirc!<chi,a6>;
d6:=Ocirc!<chi,O!(((chi`element)^-1)*Pc6)>;
//d6:=Ocirc!<mub,O!(((mub`element)^-1)*Pc6*2)>;
d6GL4:=EnhancedElementInGL4(d6);
d6GL4modN:=EnhancedElementInGL4modN(d6,N);
d6GL4^6;
EnhancedElementInGL4(d6^6);


Qz4<z4>:=CyclotomicField(4);
Zz4:=Order([1,z4]);
c4:=Embed(Zz4,O);
Pc4:=1+c4;
assert IsScalar(Pc4^4);
a4:=O!3*((chimu`element)^-1)*Pc4;
//b4:=3*Pc4/(chimu`element);
//d4:=Ocirc!<chimu,a4>;
d4:=Ocirc!<chimu,O!(((chimu`element)^-1)*Pc4)>;
d4GL4:=EnhancedElementInGL4(d4);
d4GL4modN:=EnhancedElementInGL4modN(d4,N);
d4^4;

d2GL4:=(d6GL4*d4GL4)^-1;






b:=6*Pc6/(chimu`element);
assert Norm(b) eq -6;
assert b in O;


chiM:=NormalizingElementToGL4(chi`element,O);
a6M:=UnitGroupToGL4(O!a6);
(chiM*a6M)^6;

a6rec:=EnhancedElementRecord(<chi,OmodN!a6>);
assert Order(a6rec`GL4) eq 12;
tup6:=<chi`element,a6>;
SemidirectPower(tup6,6);

brec:=EnhancedElement(<chimu,OmodN!b>);
assert Order(brec`G:4) eq 12;  //ERROR


assert Norm(a4) eq -6;
assert a4 in O;
assert Norm(b4) eq -1;
assert b4 in O;

a4rec:=EnhancedElementRecord(<chi,OmodN!a4>);
assert Order(a4rec`Gl4) eq 12; //ERROR

b4rec:=EnhancedElementRecord(<chimu,OmodN!b4>);
assert Order(b4rec`GL4) eq 8;


chirec:=EnhancedElementRecord(<chi,OmodN!(O!1)>);
assert Order(b4rec`GL4) eq 8;

//////////////////////////////////

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

ell1:=Oxfin[1,1];
c1:=Ocirc!<chimu,O!((chimu`element)^-1*(1+ell1)*3)>;
ell2:=Oxfin[2,1];
c2:=Ocirc!<chi,O!((chi`element)^-1*(1+ell2))>; 
ell3:=Oxfin[3,1];
c3:=Ocirc!<chi,O!((chi`element)^-1*(1+ell3))>;
ell4:=Oxfin[4,1];
c4:=Ocirc!<chimu,O!((chimu`element)^-1*(1+ell4)*3)>;

ell1GL4:=EnhancedElementInGL4(c1);
ell2GL4:=EnhancedElementInGL4(c2);
ell3GL4:=EnhancedElementInGL4(c3);
ell4GL4:=EnhancedElementInGL4(c4);

GL4:=Parent(ell1GL4);
assert ell1GL4^4 eq -GL4!1;
assert ell2GL4^6 eq -GL4!1;
assert ell2GL4^6 eq -GL4!1;
assert ell1GL4^4 eq -GL4!1;
assert ((ell2GL4*ell1GL4)^-2) eq GL4!1;
assert ((ell2GL4*ell4GL4)^-2) ne GL4!1;
assert ((ell3GL4*ell1GL4)^-2) eq GL4!1;
assert ((ell3GL4*ell4GL4)^-2) ne GL4!1;

CharacteristicPolynomial(((c2`element[1])`element)*(c2`element[2]));
CharacteristicPolynomial(((c3`element[1])`element)*(c3`element[2]));
//IsSimilar(ell2GL4,ell3GL4);

v6:=ell2GL4;
v4:=ell1GL4;
v2:=v6*v4;

G,Gelts:=EnhancedImageGL4(AutFull,O,N);
H:=rhocirc[10];
for H in rhocirc do
  T := CosetTable(G,H`subgroup);
  piH := CosetTableToRepresentation(G,T);
  sigma := [piH(v2),piH(v4),piH(v6)];
  Genus(sigma);
end for;


Rx<x>:=PolynomialRing(Rationals());
E<u>:=NumberField(x^2+x-1);
v:=(1-u)/2;
Rm:=Order([1,u]);

chi,map:=Embed(Rm,O);
Norm(chi);
cyc,Czeta,zeta:=IsCyclotomic(Em);


elliptic_elements_init:= [ <BxmodQx!chi,BxmodQx!(B!a)> : w in AutmuOset, x in Oxfin ];
elliptic_elements:= [ a[1]*a[2] : a in elliptic_elements_init ];


elliptic_inGL4:=[ ];

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


B<i,j,ij> := QuaternionAlgebra(Rationals(),-3,2);
O:=MaximalOrder(B);
mu := i*(2+j); assert mu^2 eq -6;
muj := mu*j/2;
//omegaO := Eltseq(K!S.2)[1] + Eltseq(K!S.2)[2]*(-1+i)/2;
//O := QuaternionOrder([1,omegaO,j,omegaO*j]);
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



Genus := function(sigma);
  //{Compute genus from permutation triple}
  d := Degree(Parent(sigma[1]));
  // Riemann-Hurwitz formula
  g2 := -2*d + &+[&+[ e[2]*(e[1]-1) : e in CycleStructure(sigma[i])] : i in [1..3]];
  assert g2 mod 2 eq 0;
  g := (g2+2) div 2;
  return g;
end function;

