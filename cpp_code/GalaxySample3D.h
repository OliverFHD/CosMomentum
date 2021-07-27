#ifndef _GalaxySample3D
#define _GalaxySample3D
#include "GalaxySample.h"


class GalaxySample3D : public GalaxySample  {

 public:

   GalaxySample3D(FlatInhomogeneousUniverseLCDM* universe, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1, BIAS_MODEL b_model);
   ~GalaxySample3D();
   
   void set_parameters_3D(double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1, BIAS_MODEL b_model);
   double set_b2_to_minimise_negative_densities_in_3D_tophat(double R_in_Mpc_over_h, double var_NL_rescale);
   double compute_variance_in_3D_tophat(double R_in_Mpc_over_h, double var_NL_rescale);
   
   int return_N_max_in_3D_tophat(double R_in_Mpc_over_h, double var_NL_rescale);
   vector<double> return_CiC_PDF_in_3D_tophat(double R_in_Mpc_over_h, double f_NL, double var_NL_rescale);
   
 private:
   
   double density; // in units with c/H_0 == 1
   double redshift;
   
};

#include "GalaxySample3D.cpp"
#endif
