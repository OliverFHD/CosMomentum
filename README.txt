
<-----> CosMomentum
<-----> README version of Dec 15 2019

<---> As of now, this code can compute
- the probability distribution function (PDF) of the matter density field,
- the cumulant generating function (CGF) of the matter density field,
- individual cumulants of the matter density field,
- all of this for Gaussian and 3 types of non-Gaussian initial conditions.
- all of this for both line-of-sight projections of the matter density field as well as the 3D density field
- all of this for PDFs of lensing convergence
- statistics of biased & stochastic tracers (e.g. galaxies)

<---> Features that will be added in the near future are
- density split statistics (cf. Gruen&Friedrich++2018, Friedrich&Gruen++2018)
- joint PDFs of pairs of cosmic random fields

<---> Installation & running:
- The notebook "compute_PDF_and_CGF.ipynb" containts examples on how to run CosMomentum.
- This includes a command that compiles the code.
- You may have to edit "cpp_code/Makefile" as suitable for your machine.

<---> Differences wrt. Friedrich et al. (2019, arXiv:1912.06621):
- This code uses the Eisenstein&Hu transfer function (incl. baryon wiggles).
- However, there is also an alternative constructor of the class Matter, that allows to read in a transfer function from a text file.

<---> How you can contribute
- Your feedback on user experience will help to improve CosMomentum!
- Please post issues or bug reports on https://github.com/OliverFHD/CosMomentum .

<---> Contributors
- the code has been designed by Oliver Friedrich, Daniel Gruen, Anik Halder, Elisabeth Krause, Tom McClintock
- we'd like to thank Cora Uhlemann for helpful discussions

<---> Acknowledgement
- Please feel free to use the code for your science and publications.
- Please cite Friedrich et al. (2019, arXiv:1912.06621) if you do so.
