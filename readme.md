# Description

This python library gathers tools and functions
related to [random vectors with a matrix representation](http://perso.quaesituri.org/florian.angeletti/Publications/Journal?filter=tags%3E[Matrix%20representation]).

In brief, this random variable are characterized by a joint probability
density function which can be written as a matrix product:
```
p(x_1,.., x_N) = tr(A^T R(x_1) ... R(x_N) ).
```

This specific library focuses on the analysis and computation
of the limit laws for the sum of such random vector.

``` 
S(X)  = x_1 + x_2 ... x_n 
```

The structure of this library follow closely the argumentation developped in the corresponding article.
It may therefore be useful for interested readers to have a look at the illustration directory.
In this directory, the illustration.py is an heavily commented python version of the relevant part 
of article whereas gen_illustration.py can be used to generate new random examples
(both generated from illustration.hyd, see [Hydra](https://github.com/FAngeletti/Hydra) ) .

For more direct use of the library, the tests avaible inside the tests directory
are quite shorter and may be more readable (at least if they were commented, it might be the case).
Note that this library relies on the sympy, numpy and matplotlib libraries. 


# Installation

The library can be installed for local user by using
```
 python setup.py install --user
```


# Script-by-script description


## Main modules 

### Connectivity.py
This module is used to compute the reachability matrix of `G(E)` and other graph theoretical quantity.
It provides for instance functions for computing the strongly connected component of `G(E)` and the
relabelling of`E` leading to block upper triangular matrix.
 
### Eigentriples.py
Classical algorithms for computing the eigenvalue and left and right eigenvectors. 
Currently, the algorithm try to use exact method for computing the eigenvalue (in order to have nice number
in the pdf file) which will of course fail for more complex system. This peculiar behavior can be turned off
by changing the default argument of the `eigentriple` function. 

### reduction.py
This module compute the reduced model for a given triple (A,E, P ( or Qs) ).

### limits.py
This module gathers the functions needed for computing the limit distributions for the generalized central limit 
theorems and law of large number from a totally irreducible model. The central limit distributions
are computed using Monte-Carlo integration. The limit laws for the (generalised) law of large number
are computed using a method based on simplicial decomposition of polytope.


### model.py
Sane entry point to the library api. This module compute the reduced model directly 
from an input triplet (A, E, Qs) where Qs are the moment matrices. Methods are then available
to compute the corresponding limit laws.


##Auxiliary modules
Modules directly used by the main module but which probably should not 
be used directly.

### polytope.py
This module implements the volume computation for general polytope in arbitrary dimension. 
The limits module use a more specific algorithm computating the volume of intersection
of polytope and hyperplane parameterizd by a polynomial. This part has by far the worst theoretical
complexity since it is factorial in the maximal length of structure chain l_max
( partially because a d dimensional hypercube is divided in d! simplices by the decomposition
algorithm currently used ). More efficient algorithm could be implemented (with a lot of free time).


### byPiece.py
This module implements piecewise polynomial and function for plotting or printing them.




## General utilities module
Modules not directly used by one of the main modules
### printing.py
Printing functions for generating readable latex file (or latex file generating readable
pdf).

### generation.py
The role of this module is to generate random model of random variable
with a matrix representation. Note that since, the interesting cases for
the limit distribution have a zero measure for most classical matrix distribution,
the library generates these models using sparse graph analysis.

### Synthesis.py
Module for generating realisations of a given m-c r.v. 
Quite slow curently.

### Histogram.py
This module implement tree-based histograms

### HistogramS.py
This module implements standard histogram
 
### utils.py 
Nameless module for miscellaneous utility functions.


