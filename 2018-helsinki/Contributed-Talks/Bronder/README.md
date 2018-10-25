The OpenCL methods to be presented at stancon 2018

To run these, download the following GPU branch and place it in a local folder

https://github.com/bstatcomp/math/tree/stancon2018

In the `_Simulations/GP.R` file add at the top:

1. The location of your cmdstan directory
2. The location of the stan math GPU library

Then you should be able to run `GP.R` to produce the results, given that your computer has OpenCL and it's drivers properly installed.

Folders:

- _Figures: The paper's figures
- _Results: Results from all of the simulations
- _Simulations: The stan and R code that run the simulations
    - _Models: The stan models
    - GP.R: Runs the simulations
    - GP_data_generator: Function to generate GP data
    - cmdstan_interface: Helper functions to connect to cmdstan from R
- report.Rmd: The report to generate Report.pdf, which contains a 4 page paper about our methods and results.


Abstract:

The Stan Math library's Hamilton Monte Carlo (HMC) sampler has computationally expensive draws while usually searching the target distribution more efficiently than alternative MCMC methods with fewer iterations. The bottleneck within draws makes Stan a prime candidate for GPU optimizations within samples. This project implements GPU optimizations for the Cholesky decomposition and it's derivative in the Stan Math library [@stanmath2015]. This work is the first known open source implementation of the Cholesky decomposition with a GPU in an HMC setting. Furthermore, the GPU kernels use OpenCL which allows the use of these methods across any brand of GPU. While results show that GPU optimizations are not optimal for small $N\times M$ matrices, large matrices can see speedups of 7.8x while retaining the same precision as models run purely on a CPU.


All the information posted here is licensed under

Code: BSD 3-clause (https://opensource.org/licenses/BSD-3-Clause)
Documentation: CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
