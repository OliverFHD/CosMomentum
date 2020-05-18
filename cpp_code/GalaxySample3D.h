#include <gsl/gsl_sf_gamma.h>


class GalaxySample3D : public GalaxySample  {

 public:

   GalaxySample3D(Matter* matter, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1);
   ~GalaxySample3D();
   
   void set_parameters_3D(double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1);
   void set_3D_bias_model_from_br_parametrisation(double b_tilde, double r, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale);
   double set_b2_to_minimise_negative_densities_in_3D_tophat(double R_in_Mpc_over_h, double var_NL_rescale);
   double compute_variance_in_3D_tophat(double R_in_Mpc_over_h, double var_NL_rescale);
   
   int return_N_max_in_3D_tophat(double R_in_Mpc_over_h, double var_NL_rescale);
   vector<double> return_CiC_PDF_in_3D_tophat(double R_in_Mpc_over_h, double f_NL, double var_NL_rescale);
   
 private:
   
   double density; // in units with c/H_0 == 1
   double redshift;
   
};
