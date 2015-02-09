README
======

Here is my MATLAB code for Tikhonov regularization with non-negativity constraint, the data processing routine used in [this paper](http://onlinelibrary.wiley.com/doi/10.1002/pssb.201451525/full). The relationship between y(t) and N_T(f) is 

```latex
y(t) = \int_0^\infty N_T(f) f \exp(-ft) df
```

where y(t) is the input data (e.g., experimental data), and N_T(f) is the spectrum to be calculated in f-domain. 


The regularization parameter can be either provided externally, or determined heuristically by L-curve criterion or Morozov discrepancy principle. Three different discretization schemes, i.e., linear, log, and Gauss-Laguerre, are available. 

Files: 

* tikregnc.m: function that solves for solutions to Tikhonov regularization with non-negativity constraint

* test_tiknc.m: the test script for tikregnc.m

* ncsolve.m: linear non-negativty constrained problem solver; a linear programming solver is used

* get_A: get discretized matrix A; linear, log, Gaussian-Laguerre quadrature discretization are available

* t_func: generates time-domain test data from a Gaussian spectrum in f-space

* gen_laguerre_rule2.m: modified Gauss-Laguerre quadrature discretizer based on [here](http://people.sc.fsu.edu/~jburkardt/m_src/gen_laguerre_rule/gen_laguerre_rule.html)

* regutools: a copy of P.C. Hansen's Regularization toolbox from [here](http://www2.imm.dtu.dk/~pcha/Regutools/); this directory should be added to MATLAB path

Other dependencies: 

* [MATLAB Optimization Toolbox](http://www.mathworks.com/products/optimization/)
