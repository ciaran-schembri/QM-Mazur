177,0
T,AlgQuatOrdRes,AlgQuatOrdResElt,0
A,AlgQuatOrdRes,2,quaternionorder,quaternionideal
A,AlgQuatOrdResElt,2,element,parent
T,AlgQuatProj,AlgQuatProjElt,0
A,AlgQuatProj,1,quaternionalgebra
A,AlgQuatProjElt,2,element,parent
S,OmodNElement,Construct an element of the OmodN whose underlying element is x in O,0,2,0,0,0,0,0,0,0,20,,0,0,AlgQuatOrdRes,,AlgQuatOrdResElt,-38,-38,-38,-38,-38
S,ElementModuloScalars,Construct an element of B^x/F^x whose underlying element is x in B,0,2,0,0,0,0,0,0,0,18,,0,0,AlgQuatProj,,AlgQuatProjElt,-38,-38,-38,-38,-38
S,eq,Decide if x equals y in OmodN,0,2,0,0,0,0,0,0,0,AlgQuatOrdResElt,,0,0,AlgQuatOrdResElt,,36,-38,-38,-38,-38,-38
S,eq,Decide if x equals y in OmodN,0,2,0,0,0,0,0,0,0,AlgQuatProjElt,,0,0,AlgQuatProjElt,,36,-38,-38,-38,-38,-38
S,eq,Decide if OmodN1 equals OmodN2,0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,AlgQuatOrdRes,,36,-38,-38,-38,-38,-38
S,eq,Decide if BxmodFx1 equals BxmodFx2,0,2,0,0,0,0,0,0,0,AlgQuatProj,,0,0,AlgQuatProj,,36,-38,-38,-38,-38,-38
S,*,compute x*y in OmodN,0,2,0,0,0,0,0,0,0,AlgQuatOrdResElt,,0,0,AlgQuatOrdResElt,,AlgQuatOrdResElt,-38,-38,-38,-38,-38
S,*,compute x*y in B^x/F^x,0,2,0,0,0,0,0,0,0,AlgQuatProjElt,,0,0,AlgQuatProjElt,,AlgQuatProjElt,-38,-38,-38,-38,-38
S,^,compute x*y in B^x/F^x,0,2,0,0,0,0,0,0,0,148,,0,0,AlgQuatProjElt,,AlgQuatProjElt,-38,-38,-38,-38,-38
S,Parent,,0,1,0,0,0,0,0,0,0,AlgQuatOrdResElt,,AlqQuatOrdRes,-38,-38,-38,-38,-38
S,Parent,,0,1,0,0,0,0,0,0,0,AlgQuatProjElt,,AlqQuatProj,-38,-38,-38,-38,-38
S,quo,,0,2,0,0,0,0,0,0,0,148,,0,0,19,,AlgQuatOrdRes,-38,-38,-38,-38,-38
S,QuaternionAlgebraModuloScalars,Create B^x/F^x,0,1,0,0,0,0,0,0,0,17,,AlgQuatProj,-38,-38,-38,-38,-38
S,IsCoercible,,0,2,0,0,0,0,0,0,0,-1,,0,0,AlgQuatOrdRes,,36,-1,-38,-38,-38,-38
S,IsCoercible,,0,2,0,0,0,0,0,0,0,-1,,0,0,AlgQuatProj,,36,-1,-38,-38,-38,-38
S,IsUnit,return whether x in O/N is a unit,0,1,0,0,0,0,0,0,0,AlgQuatOrdResElt,,36,-38,-38,-38,-38,-38
S,Set,return the set of elements O/N,0,1,0,0,0,0,0,0,0,AlgQuatOrdRes,,-50,-38,-38,-38,-38,-38
S,UnitGroup,"return (O/N)^x as a permutation group G, the second value is the isomorphism G ->(O/N)^x",0,1,0,0,0,0,0,0,0,AlgQuatOrdRes,,178,175,-38,-38,-38,-38
S,UnitGroup,"return (O/N)^x as a permutation group G, the second value is the isomorphism G ->(O/N)^x",0,2,0,0,0,0,0,0,0,148,,0,0,19,,178,175,-38,-38,-38,-38
S,ElementToAutomorphismModN,a in B^x becomes an automorphism of (O/N)^x by considering the map a |-> (x|-> a^-1xa) as long as a in N_B^x(O). We apply this to (O/N)^x as a permutation group,0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,18,,33,-38,-38,-38,-38,-38
S,ElementToAutomorphismModN,a in B^x becomes an automorphism of (O/N)^x by considering the map a |-> (x|-> a^-1xa) as long as a in N_B^x(O). We apply this to (O/N)^x as a permutation group,0,3,0,0,0,0,0,0,0,148,,0,0,19,,0,0,18,,33,-38,-38,-38,-38,-38
S,AutomorphismsModN,"Given a subset of Aut(O) input as a finite subset S of B^x, create the map theta : S -> Aut((O/N)^x)",1,0,1,83,0,AlgQuatProjElt,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,83,,175,-38,-38,-38,-38,-38
S,AutomorphismsModN,"Given a subset of Aut(O) input as a finite subset S of B^x, create the map theta : S -> Aut((O/N)^x)",1,0,1,83,0,AlgQuatProjElt,3,0,0,0,0,0,0,0,148,,0,0,19,,0,0,83,,175,-38,-38,-38,-38,-38
S,MapIsHomomorphism,Check whether the map AutmuO : C -> B^x/Q^x is an injective homomorphism,0,1,0,0,0,0,0,0,0,-1,,36,-38,-38,-38,-38,-38
S,NormalizingElementToGL4modN,O is an order over R. For an element g in N_Bx(O) the map phi_g : b |--> g^-1bg is R-linear hence [g] is an element of M_4(R) after fixing a basis this function computes [g] and also returns the R-basis of O,0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,18,,180,-38,-38,-38,-38,-38
S,NormalizingElementToGL4modN,O is an order over R. For an element g in N_Bx(O) the map phi_g : b |--> g^-1bg is R-linear hence [g] is an element of M_4(R) after fixing a basis this function computes [g] and also returns the R-basis of O,0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,AlgQuatProjElt,,180,-38,-38,-38,-38,-38
S,NormalizingElementToGL4modN,O is an order over R. For an element g in N_Bx(O) the map phi_g : b |--> g^-1bg is R-linear hence [g] is an element of M_4(R) after fixing a basis this function computes [g] and also returns the R-basis of O,0,3,0,0,0,0,0,0,0,148,,0,0,19,,0,0,18,,180,-38,-38,-38,-38,-38
S,NormalizingElementToGL4modN,O is an order over R. For an element g in N_Bx(O) the map phi_g : b |--> g^-1bg is R-linear hence [g] is an element of M_4(R) after fixing a basis this function computes [g] and also returns the R-basis of O,0,3,0,0,0,0,0,0,0,148,,0,0,19,,0,0,AlgQuatProjElt,,180,-38,-38,-38,-38,-38
S,UnitGroupModNToGL4,"O is an order over R, this returns a matrix [lambda_g] wrt to a basis which is the right regular representation lambda_g : g --> b*g where g in GL_1(O)",0,1,0,0,0,0,0,0,0,AlgQuatOrdResElt,,180,-38,-38,-38,-38,-38
S,EnhancedImagePermutation,"AutmuO is a map from a finite group C -> B^x, which is isomorphic onto the image in B^x/Q^x. We create the semidirect product of ONx by AutmuO, using AutomorphismModN as the map theta: AutmuO -> Aut(ONx)",0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,175,,-27,-38,-38,-38,-38,-38
S,EnhancedImagePermutation,"AutmuO is a map from a finite group C -> B^x, which is isomorphic onto the image in B^x/Q^x. We create the semidirect product of ONx by AutmuO, using AutomorphismModN as the map theta: AutmuO -> Aut(ONx)",0,3,0,0,0,0,0,0,0,148,,0,0,19,,0,0,-1,,-27,-38,-38,-38,-38,-38
S,EnhancedImageGL4,,0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,175,,178,-38,-38,-38,-38,-38
S,EnhancedImageGL4,,0,3,0,0,0,0,0,0,0,148,,0,0,19,,0,0,175,,178,-38,-38,-38,-38,-38
S,FixedSubspace,,0,1,0,0,0,0,0,0,0,178,,107,-38,-38,-38,-38,-38
S,HasPolarizedElementOfDegree,return an element mu of O such that mu^2 + d*disc(O) = 0 if it exists,0,2,0,0,0,0,0,0,0,148,,0,0,19,,36,20,-38,-38,-38,-38
S,DegreeOfPolarizedElement,degree of mu,0,2,0,0,0,0,0,0,0,-1,,0,0,19,,148,-38,-38,-38,-38,-38
S,IsTwisting,"(O,mu) is twisting (of degree del = -mu^2/disc(O)) if there exists chi in O and N_Bx(O) such that chi^2 = m, m|Disc(O) and mu*chi = -chi*mu",0,2,0,0,0,0,0,0,0,20,,0,0,19,,36,-38,-38,-38,-38,-38
S,Aut,,0,2,0,0,0,0,0,0,0,20,,0,0,19,,-1,-38,-38,-38,-38,-38
S,AllEnhancedSubgroups,,0,3,0,0,0,0,0,0,0,148,,0,0,20,,0,0,19,,-1,-38,-38,-38,-38,-38
S,AllEnhancedSubgroups,,0,3,0,0,0,0,0,0,0,148,,0,0,20,,0,0,17,,-1,-38,-38,-38,-38,-38
S,AllEnhancedSubgroups,,0,3,0,0,0,0,0,0,0,148,,0,0,148,,0,0,19,,-1,-38,-38,-38,-38,-38
S,AllEnhancedSubgroups,,0,3,0,0,0,0,0,0,0,148,,0,0,148,,0,0,17,,-1,-38,-38,-38,-38,-38
S,AllEnhancedSubgroups,,0,3,0,0,0,0,0,0,0,148,,0,0,148,,0,0,148,,-1,-38,-38,-38,-38,-38
S,Print,,0,1,0,0,1,0,0,0,0,AlgQuatOrdResElt,,-38,-38,-38,-38,-38,-38
S,Print,,0,1,0,0,1,0,0,0,0,AlgQuatOrdRes,,-38,-38,-38,-38,-38,-38
S,Print,,0,1,0,0,1,0,0,0,0,AlgQuatProjElt,,-38,-38,-38,-38,-38,-38
S,Print,,0,1,0,0,1,0,0,0,0,AlgQuatProj,,-38,-38,-38,-38,-38,-38
