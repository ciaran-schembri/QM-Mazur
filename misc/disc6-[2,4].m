/*rec<recformat<n: IntegerRing(), subgroup, genus, order, index, fixedsubspace, generators, split, endomorphism_representation, AutmuO_norms, ramification_data> | 
subgroup := MatrixGroup(4, IntegerRing(4)) of order 2^4
Generators:
[3 0 0 0]
[1 1 0 0]
[1 0 1 0]
[0 0 0 3]

[3 0 0 2]
[0 3 0 1]
[3 0 1 1]
[0 0 0 1]

[1 0 0 0]
[2 3 2 2]
[0 2 3 2]
[2 0 0 1]

[1 0 0 0]
[2 1 0 0]
[2 2 3 2]
[0 2 2 3],
genus := 1,
order := 16,
index := 24,
fixedsubspace := [ 2, 4 ],
generators := <<1, [1 2 2 0]>, <1, [1 2 0 2]>, <i, [1 3 3 0]>, <-3*j + k, [2 0 1 0]>>,
split := false,
endomorphism_representation := C2^2,
AutmuO_norms := { 1, 2, 3, 6 },*/

// phi is of the form X-->Y-->P^1 where Y-->P^1 is the degree 3 genus 0 Belyi map below.
ramification_data := [ Sym(24) |
(1, 2, 5, 12, 20, 11)(3, 4, 10, 17, 21, 16)(6, 7, 15, 23, 24, 14)(8, 9, 13, 22, 18, 19),
(1, 3, 8, 7)(2, 6, 14, 13)(4, 11, 18, 10)(5, 9, 16, 21)(12, 17, 22, 24)(15, 19, 20, 23),
(1, 4)(2, 7)(3, 9)(5, 13)(8, 15)(11, 19)(12, 21)(14, 22)(17, 18)(20, 24)
];

load "primitive-belyi-maps/code/primitivize.m";

X,phi:=BelyiMap(ramification_data);
X;
phi;
//////////////////////////////
//this is (X,phi)

S<x>:=PolynomialRing(Rationals());
X:=EllipticCurve(S!x^3 + 27/242*x + 189/10648,S!0);
KX<x,y> := FunctionField(X);
phi:= (864/1331*x^9 + 34992/161051*x^7 + 61236/1771561*x^6 + 472392/19487171*x^5 + 1653372/214358881*x^4 + 59049/38974342*x^3 + 11160261/25937424601*x^2 + 78121827/1141246682444*x + 182284263/50214854027536)/(x^12 + 54/121*x^10 + 54/121*x^9 + 2187/29282*x^8 + 2187/14641*x^7 + 98415/1771561*x^6 + 59049/3543122*x^5 + 38795193/3429742096*x^4 + 610173/857435524*x^3 + 129140163/207499396808*x^2 + 531441/51874849202*x + 531441/51874849202);

A,map:=MordellWeilGroup(X); //Z/2 + Z
Q:=map(A.2);  //generator
//checking ramification indices to match up with Baba-Granath's j-line we must 
//do the linear fractional transformation 0-->oo, 1-->0, oo->--16/27 which is LFT(phi,[3,2,1])

perm:=[3,2,1];
BelyiMap:=LFT(phi,perm);
a_sqfree:=SquarefreeFactorization(-6*BelyiMap);
b_sqfree:=SquarefreeFactorization(-2*(27*BelyiMap+16));


for i in [1..10] do 
  q:=i*Q;i; q;

  aq:=a_sqfree(q);
  bq:=b_sqfree(q);

  a:=TrialRepresentativeModuloSquares(aq);
  b:=TrialRepresentativeModuloSquares(bq);
  H:=QuaternionAlgebra<Rationals() | a , b >;
  IsMatrixRing(H);
end for;



MestreIsSplit:=false;
torsion:=true;
LinYang:=false;
D:=6;
exponent:=1;
bound:=50;
group:=AbelianGroup([2,4]);


//for perm in Permutations({1,2,3}) do //[3,1,2] from plugging in values looks like it has to be this one (its the only one with all 4-torsion defined over F_p^4.)
  perm:=[3,1,2];
  for i in [1..10] do 
    q:=i*Q;
    BelyiMap:=LFT(phi,perm);
    j:=BelyiMap(q);
    i;
    j;
    perm;
    MestreObstructionIsSplit(D,j : LinYang:=LinYang);
    /*if MestreIsSplit eq false or (MestreIsSplit eq true and MestreObstructionIsSplit(D,j : LinYang:=LinYang)) then  
      IgusaClebsch:=PQMIgusaClebsch(D,j : LinYang:=LinYang);
      //C:=HyperellipticCurveFromIgusaClebsch(IgusaClebsch);
      //if BaseRing(C) eq Rationals() then 
        if torsion eq true then
          tors_heur:=TorsionGroupHeuristicUpToTwist(IgusaClebsch : bound:=bound,exponent:=exponent, group:=group);
          if Type(tors_heur) eq MonStgElt then 
            tors_heur;
          else 
            invr_heur:=< PrimaryAbelianInvariants(gp) : gp in tors_heur >;
            inv:=Sprint(invr_heur);
            data:=Sprintf("%o||%o",IgusaClebsch, inv); 
            data;
          end if;
        end if; 
      //end if;
      end if;*/
  end for;
//end for;
  //GeneratePQMCurves(6 : BelyiMap:=psi, size:=10);
end for;



Rz<z>:=PolynomialRing(Rationals());
phi_num:=(32596268880384*z^9 + 10910321402112*z^7 + 1735732950336*z^6 + 1217267263872*z^5 + 387312311232*z^4 + 76079203992*z^3 + 21606265296*z^2 + 3437360388*z + 182284263);
phi_den:=(50214854027536*z^12 + 22409934855264*z^10 + 22409934855264*z^9 + 3750422981976*z^8 + 7500845963952*z^7 + 2789570813040*z^6 + 836871243912*z^5 + 568000420713*z^4 + 35734171572*z^3 + 31251919446*z^2 + 514434888*z + 514434888);


cm:=map(A.1);
phi(cm); //0



/////////////////////////////
//compute some mestre obstructions for the linear fractional transformation which .

//[1,2,3],[1,3,2],[3,2,1],[3,1,2] mestre obstruction for first 6 points.


