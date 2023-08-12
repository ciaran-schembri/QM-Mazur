//s:=GeneraTableToRecords(6,1,4 : torsioninvariants:=[2], genus:=0, endogroup:=" C2 ", sort:=true)[1];
/*
genus := 0,
fuchsindex := 8,
torsioninvariants := [ 3, 3 ],
endogroup :=  C2 ,
AutmuOnorms := { 1, 2 },
Hsplit := true,
generators :=  <<1, [2 1 1 2]>, <1, [0 0 1 1]>, <-3*j + 3*k, [0 2 2 1]>> ,
ramification_data := [
(1, 2)(3, 4)(5, 6)(7, 8),
(1, 3, 7, 6)(2, 5, 8, 4),
(1, 4)(2, 6)(3, 8)(5, 7)
]>
*/


B<i,j,k>:=QuaternionAlgebra< Rationals() | 3,-1 >;
O:=QuaternionOrder([ 1, 1/2 + 1/2*i + 1/2*j + 1/2*k, 1/2 - 1/2*i + 1/2*j - 1/2*k, 1/2 - 1/2*i - 1/2*j + 1/2*k ]);
N:=3;
Ocirc:=EnhancedSemidirectProduct(O : N:=3);

tr,mu:=HasPolarizedElementOfDegree(O,1);
mu;
assert mu^2 eq -6;
AutmuO:=Aut(O,mu);


Hgens:=[ Ocirc!<1, [2, 1, 1, 2]>, Ocirc!<1, [0, 0, 1, 1]>, Ocirc!<-3*j + 3*k, [0, 2, 2, 1]> ];
HgensGL4:=[ EnhancedElementInGL4modN(g,N) : g in Hgens ];
HgensGL4;
HGL4:=sub< GL(4,ResidueClassRing(N)) | HgensGL4 >;

FixedSubspace(HGL4);

G:=EnhancedImageGL4(AutmuO,O,N);
elliptic:=EnhancedEllipticElements(O,mu);
elliptic; 

mon:=EnhancedRamificationData(HGL4,G,O,mu);
mon;
Genus(mon);


s:=GeneraTableToRecords(6,1,3 : torsioninvariants:=[3,3], genus:=0, endogroup:=" C2 ", sort:=true)[1];
s;

/*X,phi:=BelyiMap(mon);

BelyiMap1:=LFT(phi,[3,1,2]);
BelyiMap2:=LFT(phi,[3,2,1]);


*/