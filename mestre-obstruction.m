
AttachSpec("spec");

Rx<x>:=PolynomialRing(Rationals());
//this is phi : X^*(6,P2) -> X^*(6,1)
X := Curve(ProjectiveSpace(PolynomialRing(Rationals(), 2)));
KX<x> := FunctionField(X);
//phi := KX!(32/27*x^3 - 8/9*x - 8/27);

//(Z/2)^2 torsion Belyi map
phi_init := KX!(-108*x^20 + 432*x^16 - 648*x^12 + 432*x^8 - 108*x^4)/(x^24 - 66*x^20 + 1023*x^16 + 2180*x^12 + 1023*x^8 - 66*x^4 + 1);
phi_perms:=[ LFT(phi_init,perm) : perm in Setseq(Permutations({1,2,3})) ];

phi:=phi_perms[4];

num:=Rx!Numerator(phi);
den:=Rx!Denominator(phi);

mestre_a:=(384*x^20 - 1536*x^16 + 2304*x^12 - 1536*x^8 + 384*x^4)/(x^24 + 42*x^20 + 591*x^16 + 2828*x^12 + 591*x^8 + 42*x^4 + 1);
mestre_b:=(-32*x^24 + 2112*x^20 - 32736*x^16 - 69760*x^12 - 32736*x^8 + 2112*x^4 - 32)/(x^24 + 42*x^20 + 591*x^16 + 2828*x^12 + 591*x^8 + 42*x^4 + 1);

num_a:=Rx!384*x^20 - 1536*x^16 + 2304*x^12 - 1536*x^8 + 384*x^4;
den_a:=Rx!x^24 + 42*x^20 + 591*x^16 + 2828*x^12 + 591*x^8 + 42*x^4 + 1;

num_b:=Rx!-32*x^24 + 2112*x^20 - 32736*x^16 - 69760*x^12 - 32736*x^8 + 2112*x^4 - 32;
den_b:=Rx!x^24 + 42*x^20 + 591*x^16 + 2828*x^12 + 591*x^8 + 42*x^4 + 1;

Factorization(num_a);
Factorization(den_a);
Factorization(num_b);
Factorization(den_b);


num_a1:=Rx!6;
den_a1:=Rx!(x^4 - 2*x^3 + 2*x^2 + 2*x + 1)*(x^4 + 2*x^3 + 2*x^2 - 2*x + 1);

num_b1:=Rx!-2;
den_b1:=Rx!(x^4 - 2*x^3 + 2*x^2 + 2*x + 1)*(x^4 + 2*x^3 + 2*x^2 - 2*x + 1);

for v in [20..100] do
  Discriminant(QuaternionAlgebra< Rationals() | Evaluate(num_a,v)/Evaluate(den_a,v), Evaluate(num_b,v)/Evaluate(den_b,v) >);
end for;

for v in [20..100] do
  Discriminant(QuaternionAlgebra< Rationals() | Evaluate(num_a1,v)/Evaluate(den_a1,v), Evaluate(num_b1,v)/Evaluate(den_b1,v) >);
end for;
