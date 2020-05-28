#ifndef _FlatInhomogeneousUniverseLCDM
#define _FlatInhomogeneousUniverseLCDM
#include <algorithm>    // std::min_element, std::max_element
#include <thread>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include "FlatHomogeneousUniverseLCDM.h"


class FlatInhomogeneousUniverseLCDM : public FlatHomogeneousUniverseLCDM {

 public:

  FlatInhomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max);
  FlatInhomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max, string file_for_transfer_function);
  ~FlatInhomogeneousUniverseLCDM();
  
  double return_sigma_8(){return this->cosmology.sigma_8;};
  double return_n_s(){return this->cosmology.n_s;};

  void print_Newtonian_growth_factor(string file_name);
  void return_delta_NL_of_delta_L_and_dF_ddelta_3D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values);
  void return_delta_NL_of_delta_L_and_dF_ddelta_2D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values);
  
  double Newtonian_linear_power_spectrum(double k, double e);
  double transfer_function_at(double k);
  double current_P_NL_at(double ln_k);
  double current_P_L_at(double ln_k);
  double return_linear_variance(double z, double R_in_Mpc_over_h);
  double return_non_linear_variance(double z, double R_in_Mpc_over_h);
  
  vector<double> P_L(double e);
  vector<double> P_NL(double e);

  void set_spherical_collapse_evolution_of_delta();
  void set_cylindrical_collapse_evolution_of_delta();
    
  double variance_of_matter_within_R_before_norm_was_determined(double R);
  double variance_of_matter_within_R(double R);
  double dvariance_of_matter_within_R_dR(double R);
  double variance_of_matter_within_R_NL(double R);
  double variance_of_matter_within_R_2D(double R);
  double dvariance_of_matter_within_R_dR_2D(double R);
  double variance_of_matter_within_R_NL_2D(double R);
  
  double return_D_of_z(double z);
  double return_D_of_eta(double eta);
  double return_D_prime_of_eta(double eta);
  vector<vector<double> > return_linear_growth_history(int conformal_time_steps);
  void print_growth_history(string file_name);
  
  void growth_of_DM_fluctuations_in_flat_matter_dominated_universe(double a, double *eta, double *D, double *D_prime);
  void growth_of_DM_fluctuations_in_flat_radiation_dominated_universe(double a, double *eta, double *D, double *D_prime);
  void growth_of_DM_fluctuations_in_flat_Lambda_dominated_universe(double a, double *eta, double *D, double *D_prime);
  
  vector<vector<double> > return_power_spectra(double eta, double R);
  
  vector<vector<double> > compute_phi_of_lambda_3D(double z, double R, double f_NL, double var_NL_rescale);
  void compute_phi_tilde_of_lambda_2D(double z, double R, double f_NL, vector<double> * lambda_of_delta, vector<double> * phi_of_delta, vector<double> * phi_prime_of_delta);
  vector<vector<double> > return_LOS_integrated_phi_of_lambda(double theta, double f_NL, vector<double> w_values, vector<double> n_of_w_values);
  vector<vector<vector<double> > > return_LOS_integrated_CGF_of_delta_and_kappa(double theta, double f_NL, vector<double> w_values, vector<double> n_of_w_values);
  
  vector<vector<double> > compute_PDF_3D(double z, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale);
  vector<vector<double> > compute_LOS_projected_PDF(vector<double> w_values, vector<double> n_of_w_values, double theta, double f_NL, double var_NL_rescale);
  
  int return_N_of_lambda(){return this->delta_values_for_spherical_collapse.size();};
  int return_N_of_lambda_2D(){return this->delta_values_for_cylindrical_collapse.size();};
  void return_2nd_moment_and_derivative(double R, double *variance, double *dvariance_dR);
  void return_2nd_moment_and_derivative_2D(double R, double *variance, double *dvariance_dR);
  double return_LOS_integrated_variance(double theta, vector<double> w_values, vector<double> n_of_w_values);
  double return_3D_skewness(double z, double R_in_Mpc_over_h, double f_NL);
  double return_LOS_integrated_skewness(double theta, double f_NL, vector<double> w_values, vector<double> n_of_w_values);
  
  void set_sphere_skewnesses(int PNG_modus);
  void set_sphere_skewnesses_from_eps3_powerlaw_approximation(int PNG_modus, double R_0_in_Mpc_over_h);
  void set_sphere_skewnesses_from_file(string file_name);
  
  void set_cylinder_skewnesses(int PNG_modus);
  void set_cylinder_skewnesses_from_eps3_powerlaw_approximation(int PNG_modus, double R_0_in_Mpc_over_h);
  void set_cylinder_skewnesses_from_file(string file_name);
  
  vector<double> log_top_hat_radii;
  vector<double> top_hat_sphere_variances;
  vector<double> dtop_hat_sphere_variances_dR;
  vector<double> log_top_hat_radii_for_skewnesses;
  vector<double> top_hat_sphere_skewnesses;
  vector<double> dtop_hat_sphere_skewnesses_dR;
  
  vector<double> log_top_hat_cylinder_radii;
  vector<double> top_hat_cylinder_variances;
  vector<double> dtop_hat_cylinder_variances_dR;
  vector<double> log_top_hat_cylinder_radii_for_skewnesses;
  vector<double> top_hat_cylinder_skewnesses;
  vector<double> dtop_hat_cylinder_skewnesses_dR;
  
 private:
  
  double D_initial;
  double D_prime_initial;
  double norm;
  double f_NL_rescaling_factor;
  
  vector<double> Newtonian_growth_factor_of_delta;
  vector<double> Newtonian_growth_factor_of_delta_prime;
  vector<double> Newtonian_growth_factor_second_order;
  vector<double> Newtonian_growth_factor_second_order_prime;
  
  vector<double> wave_numbers;
  vector<double> log_wave_numbers;
  vector<double> transfer_function;
  
  INITIALISATION skewness_initialisation;
  INITIALISATION skewness_initialisation_2D;

  vector<double> delta_values_for_spherical_collapse;
  vector<vector<double> > spherical_collapse_evolution_of_delta;
  vector<vector<double> > spherical_collapse_evolution_of_delta_ddelta;
  vector<vector<double> > spherical_collapse_evolution_of_delta_ddelta2;

  vector<double> delta_values_for_cylindrical_collapse;
  vector<vector<double> > cylindrical_collapse_evolution_of_delta;
  vector<vector<double> > cylindrical_collapse_evolution_of_delta_ddelta;
  vector<vector<double> > cylindrical_collapse_evolution_of_delta_ddelta2;
  
  
  void set_initial_conditions_for_growth();
  void set_wave_numbers();
  void set_transfer_function_Eisenstein_and_Hu();
  void set_transfer_function_Bond_and_Efstathiou();
  void set_transfer_function_from_file(string file);
  void set_sphere_variances();
  void set_cylinder_variances();
  
  void initialize_linear_growth_factor_of_delta();
  void initialize_up_to_second_order_growth_factor_of_delta(double D, double D_prime);  
  
   
  /*******************************************************************************************************************************************************
   * These functions and variables are for the Smith_et_al fitting formula
   *******************************************************************************************************************************************************/

  double current_k_non_linear;
  double current_n_eff;
  double current_C_sm;
  double current_scale;
   
  vector<double> current_P_NL;
  vector<double> current_P_L;
   
  double sig_sq(double R, double e);
  vector<double> c_and_n_NL(double R, double e);

  double k_NL(double k_min, double k_max, double e);

  double Delta_Q_sq(double k, double e);
  double Delta_H_sq(double k);
  double Delta_H_prime_sq(double k);
  double P_NL_at(double k, double e);
   
  void sig_sq_derivs(vector<double> (*a), vector<double> (*dadt));
  void norm_derivs(vector<double> (*a), vector<double> (*dadt));
  void c_and_n_derivs(vector<double> (*a), vector<double> (*dadt));

  double f_1, f_2, f_3, mun, nun, betan, alphan;
  double current_r_sq;
   
};

#include "FlatInhomogeneousUniverseLCDM.cpp"
#endif
