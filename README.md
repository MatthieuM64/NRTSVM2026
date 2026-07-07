# Non-reciprocal two-species Vicsek model

Codes used in the scientific publication: A. K. Dutta, S. Chatterjee, M. Mangeat, and R.Paul, <i>Stability and breakdown of chiral motion in non-reciprocal flocking</i>, <a href="https://doi.org/10.1103/yz2d-wk6g">Phys. Rev. E, accepted (2026)</a>. A preprint is available on <a href="https://arxiv.org/abs/2604.18125">arXiv</a>.</br></br>

A C++ code to compute the numerical simulations of the microscopic model is available in this repository, as well as a python code to plot the configuration of the Vicsek spins at each time. </br></br>
<b>Exportations:</b> density snapshots, particle configurations, and many observables shown in the different figures of the paper.</br>
<b>Compile:</b> g++ NRTSVM.cpp -fopenmp -lgsl -lgslcblas -lm -O3 -s -o NRTSVM.out.</br>
<b>Run:</b> ./NRTSVM.out -parameter=value.</br>
<b>List of parameters for the numerical simulations</b>: v0, eta, rho0, JAB, JBA, LX, LY, tmax, init, ran, threads (details as comments in the code).</br>
