AttachSpec("mf.spec");
probPCM := [* *];
probPQM := [* *];
probsplitPQM := [* *];
_<x> := PolynomialRing(Integers());
primeB := 100;
ell := NextPrime(80*primeB);

makemat := function(Tp, V);
  MTp := Matrix([v*Tp : v in Basis(V[1])]);
  Vmat := Matrix(Basis(V[1]));
  cols := [];
  rownum := 1;
  colnum := 1;
  while rownum le Dimension(V[1]) do
    if Vmat[rownum,colnum] eq 0 then
      colnum +:= 1;
    else
      assert Vmat[rownum,colnum] eq 1;
      Append(~cols,colnum);
      rownum +:= 1;
      colnum +:= 1;
    end if;
  end while;
  TpV := ColumnSubmatrix(MTp,cols);
  return TpV;
end function;

// for N in [2^8*3^5] do
// for N in Sort([2^i*3^j : i in [0..10], j in [0..5]]) do
// for N in Sort([c*m^2 : c in [1,2,3,6], m in [1..500] | Valuation(c*m^2,2) ne 1 and Valuation(c*m^2,3) ne 1]) do
// for N in Sort([2^10*c*m^2 : c in [1,3], m in [1..500] | m mod 2 eq 1 and Valuation(c*m^2,3) ne 1]) do
for N in [2^10*3^4] do
  // time M := NewSubspace(CuspidalSubspace(ModularSymbols(N,2,-1)));
  // print N, Factorization(N), Dimension(M);
  print "=====";
  time M := ModularSymbols(N,2,-1);
  printf "N = %o = %o, dim(M_2) = %o\n", N, Factorization(N), Dimension(M);
  if Dimension(M) eq 0 then continue; end if;
  chi0 := DirichletGroup(N)!KroneckerCharacter(1);
  for psi in TwistingCharacters(chi0,chi0) do
    // a_p(f)psi(p) = bar(a_p(f)),
    //    where bar is nontrivial involution of quadratic coefficient field QQ(sqrt(d)), d > 0
    // so psi(p) = 1 <=> a_p(f) in ZZ
    //    and |a_p(f)| <= 2*sqrt(p)
    // and psi(p) = -1 <=> a_p(f)^2 in ZZ, in fact = d*b^2
    //    and |a_p(f)|^2 <= 4*p
    print "-----";
    printf "N = %o, character = %o\n", N, psi;
    if Order(psi) eq 1 then continue; end if;
    T1 := HeckeOperator(M,1);
    T1ell := ChangeRing(SparseMatrix(T1),GF(ell));
    Vs := [[* Kernel(0*T1ell), [] *]];
    if Dimension(Vs[1][1]) eq 0 then continue; end if;
    for p in [p : p in PrimesUpTo(primeB) | N mod p ne 0 and psi(p) eq 1] do
      printf "p = %o\n", p;
      time Tp := HeckeOperator(M,p);
      time Tpspell := ChangeRing(SparseMatrix(Tp),GF(ell));
      // print p;
      newVs := [* *];
      for V in Vs do
        TpV := makemat(Tpspell, V);
        apbnd := Floor(2*Sqrt(p));
        for ap in [-apbnd..apbnd] do
          KTpap := Kernel(TpV-ap);
          if Dimension(KTpap) gt 1 then
            Append(~newVs, [* sub<V[1] | Rows(Matrix(Basis(KTpap))*Matrix(Basis(V[1])))>, V[2] cat [<p,x-ap>]*]);
          end if;
        end for;
      end for;
      Vs := newVs;
      if #Vs eq 0 then break; end if;
    end for;
    for p in [p : p in PrimesUpTo(primeB) | N mod p ne 0 and psi(p) eq -1] do
      time Tp := HeckeOperator(M,p);
      time Tpspell := ChangeRing(SparseMatrix(Tp),GF(ell));      
      // print p;
      newVs := [* *];
      for V in Vs do
        TpV := makemat(Tpspell, V);
        ap2bnd := 4*p;
        TpV2 := TpV^2;
        for ap2 in [ap2 : ap2 in [0..ap2bnd] | ap2 eq 0 or not IsSquare(ap2)] do  // trivial character, so real Hecke field
          KTpap := Kernel(TpV2-ap2);
          if Dimension(KTpap) gt 1 then
            Append(~newVs, [* sub<V[1] | Rows(Matrix(Basis(KTpap))*Matrix(Basis(V[1])))>, V[2] cat [<p,x^2-ap2>]*]);
          end if;
        end for;
      end for;
      Vs := newVs;
      if #Vs eq 0 then break; end if;
    end for;
    // eliminate CM
    ps := [p : p in PrimesUpTo(primeB) | N mod p ne 0];
    for V in Vs do
      // if Dimension(V[1]) gt 2 then continue; end if;  // probably not new
      looksCM := &or([&and([<ps[i],x> in V[2] or <ps[i],x^2> in V[2] : i in [1..#ps] | psit(ps[i]) eq -1]) : psit in TwistingCharacters(chi0,chi0) | Order(psit) ne 1]);
      if looksCM then 
        Append(~probPCM, [* N, psi, V[2], V[1] *]);
      else 
        K := SplittingField([v[2] : v in V[2]]);
        assert Degree(K) eq 2;
        d := Discriminant(K);
        r := IdentifyKroneckerCharacter(psi);
        B := QuaternionAlgebra(Rationals(), d, r);
        if IsMatrixRing(B) then
          Append(~probsplitPQM, [* N, psi, V[2], V[1] *]);
        else
          Append(~probPQM, [* N, psi, V[2], V[1], B *]);
        end if;        
      end if;
    end for;
    printf "... so far, #probPQM = %o (and #probPCM = %o and #probsplitPQM = %o)\n", 
        #probPQM, #probPCM, #probsplitPQM;
  end for;
end for;

myGcd := function(ptcnts);
  ptcnts0 := [v : v in ptcnts | v[2] ne 0];
  s1 := Gcd(ptcnts0[1][2],ptcnts0[2][2])
        *ptcnts0[1][1]^Valuation(ptcnts0[2][2],ptcnts0[1][1])
        *ptcnts0[2][1]^Valuation(ptcnts0[1][2],ptcnts0[2][1]);
  if s1 eq 1 then return 1; end if;
  return &*[p[1]^Min([Valuation(ptcnts[j][2],p[1]) : j in [1..#ptcnts] | ptcnts[j][1] ne p[1]]) : p in Factorization(s1)];  
end function;

for V in probPQM do
  N := V[1];              
  ps := [p : p in PrimesUpTo(primeB) | N mod p ne 0];
  psi := V[2];
  apchars := Sort(V[3]);
  ptcnts := [];
  for i := 1 to #ps do
    p := ps[i];
    assert apchars[i][1] eq p;
    fp := apchars[i][2];
    if psi(p) eq 1 then 
      assert Degree(fp) eq 1; 
      ap := -Evaluate(fp,0);
      Append(~ptcnts, <p,Evaluate((x^2-ap*x+p)^2,1)>);
      // print p, x^2-ap*x+p;
    else
      assert psi(p) eq -1; 
      assert Degree(fp) eq 2 and Coefficient(fp,1) eq 0; 
      ap2 := -Evaluate(fp,0);
      Append(~ptcnts, <p,Evaluate(x^4+(2*p-ap2)*x^2+p^2,1)>);
      // print p, x^4+(2*p-ap2)*x^2+p^2;
    end if;
  end for;
  print N, ptcnts, myGcd(ptcnts); // ptcnts;
end for;

/*
for V in probPQM do
  N := V[1];              
  ps := [p : p in PrimesUpTo(200) | N mod p ne 0];                                                
  chi0 := DirichletGroup(N)!KroneckerCharacter(1);
  looksCM := &or([&and([V[3][i] in [x,x^2] : i in [1..#ps] | psi(ps[i]) eq -1]) : psi in TwistingCharacters(chi0,chi0) | Order(psi) ne 1]);
  if looksCM then
    Append(~probPCM, V);
  else
    Append(~probPQMnew, V);
  end if;
end for;
*/
