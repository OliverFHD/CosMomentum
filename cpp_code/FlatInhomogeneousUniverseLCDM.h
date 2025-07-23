#ifndef _FlatInhomogeneousUniverseLCDM
#define _FlatInhomogeneousUniverseLCDM
#include <algorithm>    // std::min_element, std::max_element
#include <thread>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include "FlatHomogeneousUniverseLCDM.h"


class Pkz_grid{
  
public:
  
  Pkz_grid(vector<double> *k_values, vector<double> *z_values, vector<vector<double> > *P_L_kz_values, vector<vector<double> > *P_NL_kz_values){
    this->k_values = (*k_values);
    this->z_values = (*z_values);
    this->P_L_kz_values = (*P_L_kz_values);
    this->P_NL_kz_values = (*P_NL_kz_values);
    
    this->lnk_values = (*k_values);
    for(int i = 0; i < this->lnk_values.size(); i++){
      this->lnk_values[i] = log(this->k_values[i]);
    }
  };
  
  double return_P_L(double k, double z){
    return interpolate_neville_aitken_grid(log(k), z, &this->lnk_values, &this->z_values, &this->P_L_kz_values, constants::order_of_interpolation, constants::order_of_interpolation);
  };
  
  double return_P_NL(double k, double z){
    return interpolate_neville_aitken_grid(log(k), z, &this->lnk_values, &this->z_values, &this->P_NL_kz_values, constants::order_of_interpolation, constants::order_of_interpolation);
  };
  
private:
  
  vector<double> k_values;
  vector<double> lnk_values;
  vector<double> z_values;
  vector<vector<double> > P_NL_kz_values;
  vector<vector<double> > P_L_kz_values;
  
};


class FlatInhomogeneousUniverseLCDM : public FlatHomogeneousUniverseLCDM {

 public:

  FlatInhomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max);
  FlatInhomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max, int reduced_computation);
  FlatInhomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max, string file_for_transfer_function);
  ~FlatInhomogeneousUniverseLCDM();
  
  double return_sigma_8(){return this->cosmology.sigma_8;};
  double return_n_s(){return this->cosmology.n_s;};

  void print_Newtonian_growth_factor(string file_name);
  void return_delta_NL_of_delta_L_and_dF_ddelta_3D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values, vector<double> *delta_NL_prime_prime_values);
  void return_delta_NL_of_delta_L_and_dF_ddelta_2D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values, vector<double> *delta_NL_prime_prime_values);
  
  double Newtonian_linear_power_spectrum(double k, double e);
  double transfer_function_at(double k);
  double current_P_NL_at(double ln_k);
  double current_P_L_at(double ln_k);
  double return_linear_variance(double z, double R_in_Mpc_over_h);
  double return_non_linear_variance(double z, double R_in_Mpc_over_h, double var_NL_rescale);
  
  vector<double> P_L(double e);
  vector<double> P_NL(double e);

  void set_spherical_collapse_evolution_of_delta();
  void set_cylindrical_collapse_evolution_of_delta();
  void set_current_P_NL(double e);
    
  double variance_of_matter_within_R_before_norm_was_determined(double R);
  double variance_of_matter_within_R(double R);
  double dvariance_of_matter_within_R_dR(double R);
  double variance_of_matter_within_R_NL(double R);
  double variance_of_matter_within_R_2D(double R);
  double variance_of_matter_within_R_2D(double R, double L);
  double average_of_squared_saddle_point_within_R_2D(double R);
  double dvariance_of_matter_within_R_dR_2D(double R);
  double dvariance_of_matter_within_R_dR_2D(double R, double L);
  double d2variance_of_matter_within_R_dR2_2D(double R);
  double variance_of_matter_within_R_NL_2D(double R);
  double variance_of_matter_within_R_NL_2D(double R, double L);
  
  double return_D_of_z(double z);
  double return_D_of_eta(double eta);
  double return_D_prime_of_eta(double eta);
  vector<vector<double> > return_linear_growth_history(int conformal_time_steps);
  void print_growth_history(string file_name);
  
  void growth_of_DM_fluctuations_in_flat_matter_dominated_universe(double a, double *eta, double *D, double *D_prime);
  void growth_of_DM_fluctuations_in_flat_radiation_dominated_universe(double a, double *eta, double *D, double *D_prime);
  void growth_of_DM_fluctuations_in_flat_Lambda_dominated_universe(double a, double *eta, double *D, double *D_prime);
  
  double delta_crit(double a, double Delta);
  
  vector<vector<double> > return_power_spectra(double eta);
  
  /*
   * Computing generating functions
   * 
   */
  vector<vector<double> > compute_phi_of_lambda_3D(double z, double R, double f_NL, double var_NL_rescale);
  vector<vector<double> > compute_phi_tilde_of_lambda_2D(double e, double R, double f_NL, double var_NL_rescale);
  vector<vector<double> > return_LOS_integrated_phi_of_lambda(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values);
  vector<vector<double> > return_LOS_integrated_phi_of_lambda_incl_galaxies(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values_sources, vector<double> kernel_values_lenses);
  void return_LOS_integrated_phi_of_lambda_incl_CMB_kappa(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values, vector<vector<double> > *phi_data_delta, vector<vector<double> > *phi_data_kappa, vector<vector<double> > *phi_grid, vector<vector<double> > *dphi_dldelta_grid, vector<vector<double> > *dphi_dlkappa_grid, vector<vector<double> > *d2phi_dldelta2_grid, vector<vector<double> > *d3phi_dldelta3_grid, vector<vector<double> > *d3phi_dldelta2_dlkappa_grid, vector<vector<double> > *d3phi_dldelta_dlkappa2_grid, vector<vector<double> > *d4phi_dldelta4_grid, vector<vector<double> > *d4phi_dldelta2_dlkappa2_grid, vector<vector<double> > *d2phi_dldelta_dlkappa_grid, vector<vector<double> > *d2phi_dlkappa2_grid, vector<vector<int> > *grid_mask);
  void return_LOS_integrated_phi_of_lambda_incl_kappa(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values, vector<double> lensing_kernel_values, vector<vector<double> > *phi_data_delta, vector<vector<double> > *phi_data_kappa, vector<vector<double> > *phi_grid, vector<vector<double> > *dphi_dldelta_grid, vector<vector<double> > *dphi_dlkappa_grid, vector<vector<double> > *d2phi_dldelta2_grid, vector<vector<double> > *d2phi_dldelta_dlkappa_grid, vector<vector<double> > *d2phi_dlkappa2_grid, vector<vector<int> > *grid_mask);
  
  /*
   * Computing PDFs
   * 
   */
  vector<vector<double> > compute_PDF_3D(double z, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale);
  vector<vector<double> > compute_PDF_3D_Francis(double z, double R_in_Mpc_over_h, double f_NL, double b1, double b2, double var_NL_rescale);
  vector<vector<double> > compute_PDF_2D(double z, double R_in_Mpc_over_h, double L_in_Mpc_over_h, double f_NL, double var_NL_rescale);
  vector<vector<double> > compute_PDF_at_chosen_deltas_2D(vector<double>* chosen_deltas, double z, double R_in_Mpc_over_h, double L_in_Mpc_over_h, double f_NL, double var_NL_rescale);
  vector<vector<double> > compute_PDF_from_CGF(vector<double>* delta_values, vector<vector<double> >* CGF_data);
  vector<vector<double> > compute_PDF_from_CGF_Francis(vector<double>* delta_values, vector<vector<double> >* CGF_data);
  vector<vector<double> > compute_PDF_from_CGF_Francis_and_Gaussian_bias(vector<double>* delta_values, vector<vector<double> >* CGF_data, double b1, double b2);
  vector<vector<double> > compute_PDF_from_CGF_no_variance_terms(vector<double>* delta_values, vector<vector<double> >* CGF_data);
  vector<vector<double> > compute_LOS_projected_PDF(vector<double> w_values, vector<double> kernel_values, double theta, double f_NL, double var_NL_rescale);
  vector<vector<double> > compute_LOS_projected_PDF_incl_galaxies(vector<double> w_values, vector<double> kernel_values_sources, vector<double> kernel_values_lenses, double theta, double f_NL, double var_NL_rescale);
  vector<vector<double> > compute_LOS_projected_PDF_saddle_point(vector<double> w_values, vector<double> kernel_values, double theta, double f_NL, double var_NL_rescale, int next_to_leading);
  void compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point(double theta, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_noise_variance, vector<double> w_values, vector<double> kernel_values, vector<vector<double> > *delta_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid);
  void compute_LOS_projected_PDF_incl_kappa_saddle_point(double theta, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_noise_variance, vector<double> w_values, vector<double> kernel_values, vector<double> lensing_kernel_values, vector<vector<double> > *delta_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid);
  
  void compute_tau_coefficients(vector<double>* lambda_1st_branch, vector<double>* phi_1st_branch, vector<double>* phi_prime_1st_branch, vector<double>* phi_prime_prime_1st_branch, vector<double>* coefficients_phi_prime_of_tau, vector<double>* coefficients_dphi_prime_dtau, vector<double>* coefficients_d2phi_prime_dtau2);
  vector<vector<double> > return_LOS_joint_PDF_deltaG_kappa(double theta, double f_NL, double var_NL_rescale, double var_kappa_noise, vector<double> w_values, vector<double> kernel_values, vector<double> lensing_kernel_values);
  
  /*
   * Computing moments and cumulants
   * 
   * ISSUE: var_NL_rescale is not yet implemented in all of these functions!
   * 
   */
  void return_2nd_moment_and_derivative(double R, double *variance, double *dvariance_dR);
  //
  void return_2nd_moment_and_derivative_2D(double R, double *variance, double *dvariance_dR);
  //
  double return_3D_skewness(double z, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale);
  //
  double return_LOS_integrated_variance(double theta, vector<double> w_values, vector<double> kernel_values, double var_NL_rescale);
  //
  double return_LOS_integrated_skewness(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values);
  //
  void compute_polynomial_coefficients_from_CGF(double eta, double R, double f_NL, double var_NL_rescale, vector<double> lambda_values, vector<double> tau_values, vector<double> phi_values, vector<double> phi_prime_values, vector<double> *coeffs_phi_of_lambda, vector<double> *coeffs_phi_prime_of_lambda);
  //
  void return_LOS_integrated_polynomial_coefficients(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values, vector<double> *coeffs_phi_of_lambda, vector<double> *coeffs_phi_prime_of_lambda);
  //
  vector<double> return_LOS_integrated_C_ells(int l_max, vector<double> w_values, vector<double> kernel_values);
  //
  vector<vector<vector<double> > > return_LOS_integrated_C_ells(int l_max, vector<double> w_values, vector<vector<double> > kernel_values);
  //
  vector<vector<vector<double> > > return_LOS_integrated_3rd_moments(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<vector<double> > kernel_values);
  //
  vector<vector<double> > return_LOS_integrated_2nd_moments(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<vector<double> > kernel_values);

  
  int return_N_of_lambda(){return this->delta_values_for_spherical_collapse.size();};
  int return_N_of_lambda_2D(){return this->delta_values_for_cylindrical_collapse.size();};
  int return_Takahashi_shell_thickness(){return this->Takahashi_shell_thickness;};
  
  void set_sphere_skewnesses(int PNG_modus);
  void set_sphere_skewnesses_from_eps3_powerlaw_approximation(int PNG_modus, double R_0_in_Mpc_over_h);
  void set_sphere_skewnesses_from_file(string file_name);
  
  void compute_cylinder_skewnesses_for_unit_L_and_unit_fNL(int PNG_modus, double R_in_Mpc_over_h, double *skew, double *dskew_dR);
  void set_cylinder_skewnesses(int PNG_modus);
  void set_cylinder_skewnesses_from_eps3_powerlaw_approximation(int PNG_modus, double R_0_in_Mpc_over_h);
  void set_cylinder_skewnesses_from_file(string file_name);
  
  vector<double> log_top_hat_radii;
  vector<double> top_hat_sphere_variances;
  vector<double> dtop_hat_sphere_variances_dR;
  vector<double> d2top_hat_sphere_variances_dR2;
  vector<double> log_top_hat_radii_for_skewnesses;
  vector<double> top_hat_sphere_skewnesses;
  vector<double> dtop_hat_sphere_skewnesses_dR;
  
  vector<double> log_top_hat_cylinder_radii;
  vector<double> top_hat_cylinder_variances;
  vector<double> average_of_squared_cylinder_saddle_point;
  vector<double> dtop_hat_cylinder_variances_dR;
  vector<double> d2top_hat_cylinder_variances_dR2;
  vector<double> log_top_hat_cylinder_radii_for_skewnesses;
  vector<double> top_hat_cylinder_skewnesses;
  vector<double> dtop_hat_cylinder_skewnesses_dR;
  
  void set_power_grid(vector<double> *k_values, vector<double> *z_values, vector<vector<double> > *P_L_kz_values, vector<vector<double> > *P_NL_kz_values){
    this->power_grid = new Pkz_grid(k_values, z_values, P_L_kz_values, P_NL_kz_values);
  };
  
 private:
  
  int Takahashi_shell_thickness = 0;
  
  double D_prime_initial;
  double D_initial;
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
  
  Pkz_grid* power_grid = NULL;
  
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

  double k_NL(double k_min, double k_max, double s_min, double s_max, double e);

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
