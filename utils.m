


intrinsic getLines(file::MonStgElt) -> Any
  {turn text file into lines}
  F := Open(file, "r");
  Ls := [];
  while true do
    s := Gets(F);
    if IsEof(s) then break; end if;
    Append(~Ls,s);
  end while;
  return Ls;
end intrinsic;

intrinsic read_data(file::MonStgElt : file_type:="") -> Any
  {read the data file and make into a list with each line a record}


  if file_type eq "defining-equations" then 
    RF := recformat< n : Integers(),
    DefiningEquation,
    EndomorphismField,
    EndomorphismFieldDegree,
    EndomorphismAlgebraDiscriminantOverQ,
    TorsionSubgroupInvariants
     >;
    s := rec< RF | >;
    list:=[];
    lines:=getLines(file);
    for line in lines do
      items:=Split(line,"||");
      s`DefiningEquation:=eval(items[1]);
      s`EndomorphismField:=eval(items[2]);
      s`EndomorphismFieldDegree:=eval(items[3]);
      s`EndomorphismAlgebraDiscriminantOverQ:=eval(items[4]);
      s`TorsionSubgroupInvariants:=eval(items[5]);
      Append(~list,s);
    end for;

    return list;

  else 

    RF := recformat< n : Integers(),
    IgusaClebsch,
    TorsionHeuristicTwist
     >;
    s := rec< RF | >;
    list:=[];
    lines:=getLines(file);
    for line in lines do
      items:=Split(line,"|");
      s`IgusaClebsch:=eval(items[1]);
      s`TorsionHeuristicTwist:=eval(items[2]);
      Append(~list,s);
    end for;

    return list;
    
  end if;
end intrinsic;

intrinsic HeuristicNontrivialTorsion(file::MonStgElt : TorsionGroup:=[]) -> Any
  {return the Igusa Clebsch invariants in file which heuristically might contain TorsionGroup
  as a subgroup.}
  lines:=read_data(file);
  data:=[];
  if TorsionGroup eq [] then
    for s in lines do
      Append(~data,s);
    end for;
  else
    A:=AbelianGroup(TorsionGroup);
    for s in lines do
      possible_tors:=[ AbelianGroup(abinv) : abinv in s`TorsionHeuristicTwist | abinv ne [] ];
      if true in [ IsSubgroup(A,B) : B in possible_tors ] then
        Append(~data,s);
      end if;
    end for;
  end if;
  return data;
end intrinsic;


intrinsic GlobalPQMFromHeuristics(file::MonStgElt : TorsionGroup:=[], bound:=1000, ComputeEndomorphisms:=false, writetofile:="",TwistSearch:=false) -> Any 
  {}
  Rx<x>:=PolynomialRing(Rationals());
  prec := 1000;
  F := RationalsExtra(prec);
  list:=HeuristicNontrivialTorsion(file: TorsionGroup:=TorsionGroup);

  for s in list do
    try 
      IgusaClebsch:=ChangeUniverse(s`IgusaClebsch,Rationals());
      inv:=s`TorsionHeuristicTwist;
      C:=HyperellipticCurveFromIgusaClebsch(IgusaClebsch);
      X:=ReducedWamelenModel(C);
      if TwistSearch eq true then 
        Xtors:=NaiveTorsionSearchTwist(X,TorsionGroup: bound:=bound);
      else 
        Xtors:=X;
      end if;

      if Type(Xtors) eq CrvHyp then
        print "=====================";
        //Factorization(Discriminant(Xtors));
        tors_gp:=PrimaryAbelianInvariants(TorsionSubgroup(Jacobian(Xtors)));
        
        XF := ChangeRing(Xtors,F);
        if ComputeEndomorphisms eq true then 
          _,B:=HeuristicEndomorphismAlgebra(XF : CC:=true);
          tr,D:=IsQuaternionAlgebra(B);
          Discriminant(D);
        end if;

    
        _,E:=HeuristicEndomorphismAlgebra(XF : Geometric:=false);
        _,E:=IsNumberField(E);
        L:=HeuristicEndomorphismFieldOfDefinition(Xtors);
        //assert Degree(L) eq 4;
        pols:=Sprint(Xtors,"Magma");
        fld:= Sprint(DefiningPolynomial(L),"Magma");
        data:=Sprintf("%o||%o||%o||%o||%o",pols, fld,Degree(L),Discriminant(E),tors_gp); data;
        if writetofile ne "" then 
          PrintFile(writetofile,data);
        end if;
      end if;
    catch e
     ;
    end try;
  end for;
  return "";
end intrinsic;

intrinsic IsNumberField(A::AlgAss) -> BoolElt, FldNum
  {Take the associative algebra over Q and make a number field}
  if BaseRing(A) eq Rationals() and IsCommutative(A) then 
    B:=Basis(A); 
    pols:=[ MinimalPolynomial(b) : b in B ];
    F:=NumberField(pols);
    return true,F;
  else 
   return false;
  end if;
end intrinsic;


intrinsic SquarefreeFactorization(x::FldRatElt) -> FldRatElt
  {squarefree factorization of a rational number}
  numx:=Numerator(x);
  denx:=Denominator(x);

  numa,c:=SquarefreeFactorization(numx);
  dena,d:=SquarefreeFactorization(denx);

  return numa/dena;
end intrinsic;


intrinsic RepresentativeModuloSquares(x::FldRatElt) -> RngIntElt
  {x = q^2*a where a is a squarefree integer, return a}

  numx:=Numerator(x);
  denx:=Denominator(x);

  numa,c:=SquarefreeFactorization(numx);
  dena,d:=SquarefreeFactorization(denx);
  
  a:=numa*dena;
  assert IsSquarefree(a);
  assert IsSquare(x/a);
  return numa*dena;
end intrinsic;


intrinsic TrialRepresentativeModuloSquares(x::FldRatElt : divisionbound:=1000000000) -> RngIntElt
  {x = q^2*a where a is a squarefree integer, return a}

  numx:=Numerator(x);
  denx:=Denominator(x);

  trinum1,trinum2:=TrialDivision(numx,divisionbound);
  trinum2:=trinum2 cat [1];
  assert (&*[ a[1]^a[2] : a in trinum1 ])*trinum2[1] eq Abs(numx);

  triden1,triden2:=TrialDivision(denx,divisionbound);
  triden2:=triden2 cat [1];
  assert (&*[ a[1]^a[2] : a in triden1 ])*triden2[1] eq Abs(denx);

  newnum:=(&*[ a[1]^(a[2] mod 2) : a in trinum1 ])*trinum2[1];
  newden:=(&*[ a[1]^(a[2] mod 2) : a in triden1 ])*triden2[1];

  newx:=Sign(x)*newnum*newden;

  return newx;
end intrinsic;




intrinsic SquarefreeFactorization(phi::FldFunFracSchElt[CrvEll[FldRat]]) -> FldFunFracSchElt[CrvEll[FldRat]]
  {If phi(x) = f(x)/g(x) and f = c^2*f0, g = d^2*g0 with f0, g0 irreducible; return f0/g0.}
  
  KX<x,y>:=Parent(phi);
  phi_num:=KX!Numerator(phi);
  phi_den:=KX!Denominator(phi);

  Rz<z>:=PolynomialRing(Rationals());
  phi_numz:=Rz!Evaluate(phi_num,[z]);
  phi_denz:=Rz!Evaluate(phi_den,[z]);

  fac_numz:=Factorization(phi_numz);
  fac_denz:=Factorization(phi_denz);

  fac_l1:=Factorization(LeadingCoefficient(phi_numz));
  fac_l2:=Factorization(LeadingCoefficient(phi_denz));

  newnumz:=(&*[ a[1]^(a[2] mod 2) : a in fac_l1 ])*(&*[ a[1]^(a[2] mod 2) : a in fac_numz ]);
  newdenz:=(&*[ a[1]^(a[2] mod 2) : a in fac_l2 ])*(&*[ a[1]^(a[2] mod 2) : a in fac_denz ]);
  
  phi_sqfree:=(KX!Evaluate(newnumz,x))/(KX!Evaluate(newdenz,x));

  return phi_sqfree;
end intrinsic;






