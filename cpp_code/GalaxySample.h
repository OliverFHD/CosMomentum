#include <gsl/gsl_sf_gamma.h>


class GalaxySample {

 public:

  GalaxySample(Matter* matter, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1);
  ~GalaxySample();
  
  Universe* universe;
  Matter* matter;
  
  void change_parameters(double density_in_Mpc_over_h_cubed, double z, double b1, double b2, double a0, double a1);
  double set_b2_to_minimise_negative_densities(double z, double R_in_Mpc_over_h);
  
  int return_N_max(double z, double R_in_Mpc_over_h);
  int return_N_max_and_variance(double z, double R_in_Mpc_over_h, double* variance);
  
  vector<double> return_CiC_PDF(double z, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale);
  
 private:
   
   double density; // int units with c/H_0 == 1
   double redshift;
   double linear_bias;
   double quadratic_bias;
   double alpha_0; // shot-noise parameters a la https://arxiv.org/abs/1710.05162
   double alpha_1; // shot-noise parameters a la https://arxiv.org/abs/1710.05162
   // --> ISSUE: this pretends that quadratic_bias, alpha_0 and alpha_1 can be the same on all scales. That is 
   //            however not possible. Scale dependence should be made explicit in the future!
   
  cosmological_model cosmology;
   
   int error_flag_negative_density = 0;
   
   double return_P_of_N_given_delta(int N, double V, double delta, double variance); // variance only used when quadratic_bias != 0
   
};
