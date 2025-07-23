
<-----> CosMomentum
<-----> README version of July 23 2025


*******************************************
IF YOU ARE LOOKING FOR CODE TO SIMULATE 
OVERLAPPING QUBITS / RESULTS FROM
"Holographic phenomenology via overlapping degrees of freedom" 
https://arxiv.org/abs/2402.11016
PLEASE CHANGE TO THIS REPOSITORY:

https://github.com/OliverFHD/GPUniverse

(There was a wrong link in the arXiv submission...)
*******************************************

<---> As of now, this code can compute
- the probability distribution function (PDF) of the matter density field,
- the cumulant generating function (CGF) of the matter density field,
- individual cumulants of the matter density field,
- all of this for Gaussian and 3 types of non-Gaussian initial conditions.
- all of this for both line-of-sight projections of the matter density field as well as the 3D density field
- all of this for PDFs of lensing convergence
- statistics of biased & stochastic tracers (e.g. galaxies)
- joint PDF of galaxy density and lensing convergence


<---> Installation & running:
- The notebook "compute_1D_PDFs_and_CGFs.ipynb" containts examples on how to run CosMomentum for single-field statistics.
- The notebook "compute_2D_PDF.ipynb" containts examples on how to run CosMomentum for the joint PDF of galaxy density and lensing convergence.
- This includes a command that compiles the code.
- You may have to edit "cpp_code/Makefile" as suitable for your machine.

<---> Differences wrt. Friedrich et al. (2019, arXiv:1912.06621):
- This code uses the Eisenstein&Hu transfer function (incl. baryon wiggles).
- However, there is also an alternative constructor of the class Matter, that allows to read in a transfer function from a text file.

<---> How you can contribute
- Your feedback on user experience will help to improve CosMomentum!
- Please post issues or bug reports on https://github.com/OliverFHD/CosMomentum .

<---> Contributors
- the code has been designed by Oliver Friedrich, Daniel Gruen, Anik Halder, Lina Castiblanco, Cora Uhlemann, Elisabeth Krause, Tom McClintock

<---> Acknowledgement
- Please feel free to use the code for your science and publications.
- Please cite Friedrich et al. (2019, arXiv:1912.06621) if you do so.
