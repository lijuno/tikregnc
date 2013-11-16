README
======

Here are my MATLAB sample codes for Tikhonov regularization with non-negativity constraint. The relationship between y(t) and N_T(f) is 

y(t) = \int_0^\infty N_T(f) f \exp(-ft) df

where y(t) is the input data (e.g., experimental data), and N_T(f) is the spectrum to be calculated in f-domain. 

Files: 

* tikregnc.m: function that solves for solutions to Tikhonov regularization with non-negativity constraint

* test_tiknc.m: the test script for tikregnc.m

* ncsolve.m: linear non-negativty constrained problem solver

* get_A: get discretized matrix A

* t_func: generates time-domain test data from a Gaussian spectrum in f-space

* gen_laguerre_rule2.m: modified Gauss-Laguerre quadrature discretizer based on [here](http://people.sc.fsu.edu/~jburkardt/m_src/gen_laguerre_rule/gen_laguerre_rule.html)

* regutools: a copy of P.C. Hansen's Regularization toolbox from [here](http://www2.imm.dtu.dk/~pcha/Regutools/); this directory should be added to MATLAB path

Other dependencies: 

* [MATLAB Optimization Toolbox](http://www.mathworks.com/products/optimization/)