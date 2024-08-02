Below is auxilliary code for the paper: 

"RATIONAL TORSION POINTS ON ABELIAN SURFACES WITH QUATERNIONIC MULTIPLICATION" 
https://arxiv.org/pdf/2308.15193

authors: Jef Laga, Ciaran Schembri, Ari Shnidman, John Voight

To verify the calculations in remark 2.3.3, either copy-paste the contents of the following file or load into Magma with:

> load "quaternionfixedpoints.m";

The code used for the modular form computations in proposition 5.1.7 is in the file "qmsearch.m"

To load the defining equation f of genus 2 curve $C_t$ in proposition 5.3.9, open Magma and type:

> load "prop5-3-9-Z2Z2-defining-equations.m"; <br />
> f;

To verify the computations in the proof of proposition 5.3.9, one can copy-paste into Magma or load the file:
load "prop5-3-9-Z2Z2.m";

The principally polarized PQM examples in Table 2 can be found in the file "table2-principally-polarized-examples.m"


