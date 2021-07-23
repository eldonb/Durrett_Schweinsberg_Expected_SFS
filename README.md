# Durrett_Schweinsberg_Expected_SFS

The C code computes exact expected branch length spectrum for 
Example 2.4 in Durrett and Schweinsberg (2005). 

You will need the GSL library. The code should be error-free, and you can compile it with e.g.

 gcc -Wall -O3 -mtune=corei7 -march=native -DNDEBUG ebi_eb_logits_DS_Example24.c -lm -lgsl -lgslcblas

and  run it as 
./a.out (sample size)  (c parameter)
