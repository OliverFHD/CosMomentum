#include <thread>
#include <random>
#include <stdio.h>
#include <stdlib.h>


class Matter {

 public:

  Matter(Universe* uni);
  Matter(Universe* uni, string file_for_transfer_function);
  ~Matter();
  
  Universe* universe;

  void print_Newtonian_growth_factor(string file_name);
  void return_delta_NL_of_delta_L_and_dF_ddelta_3D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values);
  void return_delta_NL_of_delta_L_and_dF_ddelta_2D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values);
  
  double Newtonian_linear_power_spectrum(double k, double e);
  double transfer_function_at(double k);
  double current_P_NL_at(double ln_k);
  double current_P_L_at(double ln_k);
  double P_L_today_at(double ln_k);
  double return_linear_variance(double z, double R_in_Mpc_over_h);
  double return_non_linear_variance(double z, double R_in_Mpc_over_h);
  void return_2D_non_linear_variance(double *var_deltaLOS, double *var_GaussLOS, double *var_tophatLOS, double z, double R_in_Mpc_over_h, double L_in_Mpc_over_h);
  
  vector<double> P_L(double e);
  vector<double> P_NL(double e);
  vector<double> return_wave_numbers();
  
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

  void print_P_NL(double w, string output_file);
  void set_spherical_collapse_evolution_of_delta(double z_min, double z_max, int n_time);
  void set_cylindrical_collapse_evolution_of_delta(double z_min, double z_max, int n_time);
    
  double variance_of_matter_within_R_before_norm_was_determined(double R);
  double variance_of_matter_within_R(double R);
  double dvariance_of_matter_within_R_dR(double R);
  double variance_of_matter_within_R_NL(double R);
  double variance_of_matter_within_R_2D(double R);
  double variance_of_matter_within_R_2D(double R, double L);
  double variance_of_matter_within_R_2D_GaussianLOS(double R, double L);
  double dvariance_of_matter_within_R_dR_2D(double R);
  double dvariance_of_matter_within_R_dR_2D(double R, double L);
  double dvariance_of_matter_within_R_dR_2D_GaussianLOS(double R, double L);
  double variance_of_matter_within_R_NL_2D(double R);
  double variance_of_matter_within_R_NL_2D(double R, double L);
  double variance_of_matter_within_R_NL_2D_GaussianLOS(double R, double L);
  
  ///// Needed for PNG calculation:
  double skewness_of_matter_within_R(double R, double alpha_1, double alpha_2, double alpha_3);
  double dskewness_of_matter_within_R_dR(double R, double alpha_1, double alpha_2, double alpha_3);
  /////
  
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
  vector<vector<double> > compute_phi_of_lambda_2D(double z, double R, double L, double f_NL, double var_NL_rescale, int LOS_modus);
  void compute_phi_tilde_of_lambda_2D(double z, double R, double f_NL, vector<double> * lambda_of_delta, vector<double> * phi_of_delta, vector<double> * phi_prime_of_delta);
  vector<vector<double> > compute_phi_of_lambda_3D_EdS(double z, double R, double f_NL, double var_NL_rescale);
  vector<vector<double> > return_LOS_integrated_phi_of_lambda(double theta, double f_NL, vector<double> z_values, vector<double> n_of_z_values);
  
  vector<vector<double> > compute_PDF_3D(double z, double R, double f_NL, double var_NL_rescale);
  vector<vector<double> > compute_LOS_projected_PDF(vector<double> z_values, vector<double> n_of_z_values, double theta, double f_NL, double var_NL_rescale);
  
  int return_N_of_lambda(){return this->delta_values_for_spherical_collapse.size();};
  int return_N_of_lambda_2D(){return this->delta_values_for_cylindrical_collapse.size();};
  void return_2nd_moment_and_derivative(double R, double *variance, double *dvariance_dR);
  void return_2nd_moment_and_derivative_2D(double R, double *variance, double *dvariance_dR);
  double return_LOS_integrated_variance(double theta, vector<double> z_values, vector<double> n_of_z_values);
  
  void set_sphere_skewnesses(int PNG_modus);
  void set_sphere_skewnesses_from_eps3_powerlaw_approximation(int PNG_modus, double R_0_in_Mpc_over_h);
  void set_sphere_skewnesses_from_file(string file_name);
  
  void set_cylinder_skewnesses(int PNG_modus);
  void set_cylinder_skewnesses_from_eps3_powerlaw_approximation(int PNG_modus, double R_0_in_Mpc_over_h);
  void set_cylinder_skewnesses_from_file(string file_name);
  
 private:
   
  int number_of_entries_Newton;
  int ell_max;

  double a_initial;
  double a_final;
  double eta_initial;
  double eta_final;
  double D_initial;
  double D_prime_initial;
  double norm;
  double f_NL_rescaling_factor;
  double top_hat_radius;
  double second_top_hat_radius;
  
  INITIALISATION skewness_initialisation;
  INITIALISATION skewness_initialisation_2D;
      
  vector<double> eta_Newton;   
  vector<double> eta_NL;   
  vector<double> Newtonian_growth_factor_of_delta;
  vector<double> Newtonian_growth_factor_of_delta_prime;
  vector<double> Newtonian_growth_factor_of_delta_prime_prime;
  vector<double> Newtonian_growth_factor_of_delta_prime_prime_prime;
  vector<double> Newtonian_growth_factor_second_order;
  vector<double> Newtonian_growth_factor_second_order_prime;
  
  vector<double> wave_numbers;
  vector<double> log_wave_numbers;
  vector<double> transfer_function;
  vector<double> P_L_today;
  vector<double> P_NL_today;

  vector<double> delta_values_for_spherical_collapse;
  vector<vector<double> > spherical_collapse_evolution_of_delta;
  vector<vector<double> > spherical_collapse_evolution_of_delta_ddelta;
  vector<vector<double> > spherical_collapse_evolution_of_delta_ddelta2;
  vector<double> eta_NL_for_spherical_collapse;

  vector<double> delta_values_for_cylindrical_collapse;
  vector<vector<double> > cylindrical_collapse_evolution_of_delta;
  vector<vector<double> > cylindrical_collapse_evolution_of_delta_ddelta;
  vector<vector<double> > cylindrical_collapse_evolution_of_delta_ddelta2;
  vector<double> eta_NL_for_cylindrical_collapse;
  
  cosmological_model cosmology;
  
  void set_initial_conditions();
  void set_wave_numbers();
  void set_transfer_function_Eisenstein_and_Hu();
  void set_transfer_function_Bond_and_Efstathiou();
  void set_transfer_function_from_file(string file);
  void set_P_today();
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







void return_EdS_spherical_collapse(double D, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values){
  
  double delta_min = -10.0;
  double delta_max = 1.674;
  double ddelta = 0.001;
  double delta = delta_min-ddelta;
  
  int n = 1+int((delta_max - delta_min)/ddelta);
  
  
  delta_L_values->resize(n, 0.0);
  delta_NL_values->resize(n, 0.0);
  delta_NL_prime_values->resize(n, 0.0);
  
  for(int i = 0; i < n; i++){
    delta += ddelta;
    (*delta_L_values)[i] = delta;
  }
  
  fstream input;
  input.open("Data/pretabulated_data_for_EdS/EdS_spherical_collapse.dat");
  
  string headline;
  double dummy;
  vector<double> dummy_delta_L(0, 0.0);
  vector<double> dummy_delta_NL(0, 0.0);
  vector<double> dummy_delta_NL_prime(0, 0.0);
  
  std::getline(input, headline);
  
  while(input.good()){
    input >> dummy;
    if(input.good()){
      dummy_delta_L.push_back(dummy);
      input >> dummy; dummy_delta_NL.push_back(dummy);
      input >> dummy; dummy_delta_NL_prime.push_back(dummy);
    }
  }
  
  
  for(int i = 0; i < n; i++){
    delta = D*(*delta_L_values)[i];
    (*delta_NL_values)[i] = interpolate_neville_aitken(delta, &dummy_delta_L, &dummy_delta_NL, constants::order_of_interpolation);
    (*delta_NL_prime_values)[i] = D*interpolate_neville_aitken(delta, &dummy_delta_L, &dummy_delta_NL_prime, constants::order_of_interpolation);
  }
  
  input.close();
}

