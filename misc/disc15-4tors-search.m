
Rx<x>:=PolynomialRing(Rationals());
//this is phi : X^*(6,P2) -> X^*(6,1)
//X := Curve(ProjectiveSpace(PolynomialRing(Rationals(), 2)));
//KX<x> := FunctionField(X);
//phi1Elkies := Rx!(-1*(16*x^6-24*x^4+9*x^2-1))-1;
//phi1LY:=(-1/27)*phi1Elkies;
//BelyiMap:=phi1LY;


phiX2Elkies:=Rx!x*(x-3)^2/4;
c:=1;
phiX4Elkies:=Rx!Evaluate(phiX2Elkies,c*x^2);
phiX4LY:=(-1/27)*phiX4Elkies;
BelyiMap:=phiX4LY;

D:=15;
LinYang:=true;


height_init:=1;
size:=1000;
torsion_jsQ:=[];
igusa_clebsch_list:=[];
size:=size;
height:=height_init;
while #torsion_jsQ lt size do
  height:=height+1;
  extra_a:=[ height/b : b in [-height..height] | b ne 0 and GCD(height,b) eq 1 ];
  extra_b:=[ a/height : a in [-height..height] | a ne 0 and GCD(height,a) eq 1 ];
  extra := extra_a cat extra_b;
  for q in extra do

  	try

	    j:=Evaluate(BelyiMap,q);
	    IgusaClebsch:=PQMIgusaClebsch(D,j : LinYang:=LinYang);
	 
	    Append(~torsion_jsQ,j);
	    if MestreObstructionIsSplit(D,j) then 
	      j;   Append(~igusa_clebsch_list,IgusaClebsch);

	    end if;
	catch e 
	  ;
	end try;

  end for;
end while;



sprint:=[ #Sprint(I) : I in igusa_clebsch_list ];
ParallelSort(~sprint,~igusa_clebsch_list);
/*for i in [1..#igusa_clebsch_list] do 
  i;
  #Sprint(HyperellipticCurveFromIgusaClebsch(igusa_clebsch_list[i]));
end for;*/


for IgusaClebsch in igusa_clebsch_list do 
  TorsionGroupHeuristicUpToTwist(IgusaClebsch : group:=AbelianGroup([4]), bound:=500, exponent:=2);
	//C:=HyperellipticCurveFromIgusaClebsch(IgusaClebsch);
	//X:=ReducedWamelenModel(C);
  //PrimaryAbelianInvariants(TorsionSubgroup(Jacobian(X)));
end for;

for IgusaClebsch in igusa_clebsch_list do 

	C:=HyperellipticCurveFromIgusaClebsch(IgusaClebsch);
 if BaseRing(C) eq Rationals() then 
    if BaseRing(C) eq Rationals() then 
      X:=ReducedWamelenModel(C);
    else 
      X:=C;
    end if;
    X;

    prec := 1000;
    if BaseRing(X) eq Rationals() then 
    	F:=RationalsExtra(prec);
        XF:=ChangeRing(X,F);
    else 
        F := BaseNumberFieldExtra(DefiningPolynomial(BaseRing(X)),prec);
        tr,map:=IsIsomorphic(BaseRing(C),F);

        f:=HyperellipticPolynomials(X);
        coefs:= [map(a) : a in Coefficients(f)];
        XF:=HyperellipticCurve(coefs);
    end if;

    K:=HeuristicEndomorphismFieldOfDefinition(XF); K;
    XK:=ChangeRing(X,K);
    TorsionGroupHeuristicUpToTwist(XK : bound:=1000);

    _,B:=HeuristicEndomorphismAlgebra(XF : CC:=true);
    B;
    tr,A:=IsQuaternionAlgebra(B);
    tr;
    if tr then Discriminant(A); end if;
  end if;
end for;

