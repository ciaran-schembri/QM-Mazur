177,0
T,AlgQuatOrdRes,AlgQuatOrdResElt,0
A,AlgQuatOrdRes,2,quaternionorder,quaternionideal
A,AlgQuatOrdResElt,2,element,parent
S,OmodNElement,Construct an element of the OmodN whose underlying element is x in O,0,2,0,0,0,0,0,0,0,20,,0,0,AlgQuatOrdRes,,AlgQuatOrdResElt,-38,-38,-38,-38,-38
S,eq,Decide if x equals y in OmodN,0,2,0,0,0,0,0,0,0,AlgQuatOrdResElt,,0,0,AlgQuatOrdResElt,,36,-38,-38,-38,-38,-38
S,eq,Decide if OmodN1 equals OmodN2,0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,AlgQuatOrdRes,,36,-38,-38,-38,-38,-38
S,*,compute x*y in OmodN,0,2,0,0,0,0,0,0,0,AlgQuatOrdResElt,,0,0,AlgQuatOrdResElt,,36,-38,-38,-38,-38,-38
S,Parent,,0,1,0,0,0,0,0,0,0,AlgQuatOrdResElt,,AlqQuatOrdRes,-38,-38,-38,-38,-38
S,quo,,0,2,0,0,0,0,0,0,0,148,,0,0,19,,AlgQuatOrdRes,-38,-38,-38,-38,-38
S,IsCoercible,,0,2,0,0,0,0,0,0,0,-1,,0,0,AlgQuatOrdRes,,36,-1,-38,-38,-38,-38
S,IsUnit,return whether x in O/N is a unit,0,1,0,0,0,0,0,0,0,AlgQuatOrdResElt,,36,-38,-38,-38,-38,-38
S,Set,return the set of elements O/N,0,1,0,0,0,0,0,0,0,AlgQuatOrdRes,,-50,-38,-38,-38,-38,-38
S,UnitGroup,"return (O/N)^x as a permutation group G, the second value is the isomorphism G ->(O/N)^x",0,1,0,0,0,0,0,0,0,AlgQuatOrdRes,,224,175,-38,-38,-38,-38
S,UnitGroupQuotient,"return (O/N)^x as a permutation group G, the second value is the isomorphism G ->(O/N)^x",0,2,0,0,0,0,0,0,0,148,,0,0,19,,224,175,-38,-38,-38,-38
S,ElementToAutomorphismModN,a in B^x becomes an automorphism of (O/N)^x by considering the map a |-> (x|-> a^-1xa) as long as a in N_B^x(O). We apply this to (O/N)^x as a permutation group,0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,18,,33,-38,-38,-38,-38,-38
S,ElementToAutomorphismModN,a in B^x becomes an automorphism of (O/N)^x by considering the map a |-> (x|-> a^-1xa) as long as a in N_B^x(O). We apply this to (O/N)^x as a permutation group,0,3,0,0,0,0,0,0,0,148,,0,0,19,,0,0,18,,33,-38,-38,-38,-38,-38
S,AutomorphismModN,"Given a subgroup of Aut(O) input as a finite subset A of B^x, create the map theta : A -> Aut((O/N)^x)",0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,-1,,175,-38,-38,-38,-38,-38
S,EnhancedImage,"AutmuO is a map from a finite group C -> B^x, which is isomorphic onto the image in B^x/Q^x. We create the semidirect product of ONx by AutmuO, using AutomorphismModN as the map theta: AutmuO -> Aut(ONx)",0,2,0,0,0,0,0,0,0,AlgQuatOrdRes,,0,0,-1,,-27,-38,-38,-38,-38,-38
S,EnhancedImage,"AutmuO is a map from a finite group C -> B^x, which is isomorphic onto the image in B^x/Q^x. We create the semidirect product of ONx by AutmuO, using AutomorphismModN as the map theta: AutmuO -> Aut(ONx)",0,3,0,0,0,0,0,0,0,148,,0,0,19,,0,0,-1,,-27,-38,-38,-38,-38,-38
S,Print,,0,1,0,0,1,0,0,0,0,AlgQuatOrdResElt,,-38,-38,-38,-38,-38,-38
S,Print,,0,1,0,0,1,0,0,0,0,AlgQuatOrdRes,,-38,-38,-38,-38,-38,-38
