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

ramification_data := [ Sym(24) |
(1, 2, 5, 12, 20, 11)(3, 4, 10, 17, 21, 16)(6, 7, 15, 23, 24, 14)(8, 9, 13, 22, 18, 19),
(1, 3, 8, 7)(2, 6, 14, 13)(4, 11, 18, 10)(5, 9, 16, 21)(12, 17, 22, 24)(15, 19, 20, 23),
(1, 4)(2, 7)(3, 9)(5, 13)(8, 15)(11, 19)(12, 21)(14, 22)(17, 18)(20, 24)
];

GG:=sub< Sym(24) | ramification_data >;
IsPrimitive(GG);
X,phi:=BelyiMap(ramification_data);
X;
phi;



