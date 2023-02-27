SetColumns(0);
Rx<x>:=PolynomialRing(Rationals());
BelyiMap:=-1/27*(x*(x-3)^2/4);

GeneratePQMCurves(15 : size:=10000, endomorphisms:=false, BelyiMap:=BelyiMap);

/*
GlobalPQMFromHeuristics("PQMdata/data/X*(15,2)-igusaclebsch1.m":
	writetofile:="PQMdata/data/X*(15,2)-defining-equations1.m");
*/