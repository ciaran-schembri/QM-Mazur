

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


////////////////////////////////////

B:=QuaternionAlgebra(6);
BxmodQx:=QuaternionAlgebraModuloScalars(B);
O:=MaximalOrder(B);
Ocirc:=EnhancedSemidirectProduct(O);
N:=3;

del:=1;
tr,mu:=HasPolarizedElementOfDegree(O,del);

rhocirc:=AllEnhancedSubgroups(O,mu,N: minimal:=true, verbose:=true, PQMtorsion:=true);


