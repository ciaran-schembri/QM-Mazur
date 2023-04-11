

B:=QuaternionAlgebra(14);
O:=MaximalOrder(B);
N:=4;

SetProfile(true);
for del in Divisors(Discriminant(B)) do
  if HasPolarizedElementOfDegree(O,del) then 
    tr,mu:=HasPolarizedElementOfDegree(O,del);
    rhocirc:=AllEnhancedSubgroups(O,mu,N: minimal:=true, verbose:=true, onlypotential:=true);
    G:=ProfileGraph();
    ProfilePrintByTotalTime(G : Max:=10);
  end if;
 end for;




