
#include "constants.h"
#include "lognormal_tools.h"
#include "io_utils.h"
#include "error_handling.h"
#include "interpolators.h"
#include "LA_utils.h"
#include "FlatInhomogeneousUniverseLCDM.h"
#include "GalaxySample3D.h"
#include "ProjectedGalaxySample.h"
//#include "/Library/Frameworks/Python.framework/Versions/3.7/include/python3.7m/Python.h"


class Global_Universes {

  public:
    Global_Universes(){}
    
    cosmological_model cosmo;
    
    vector<vector<double> > delta_grid;
    vector<vector<double> > kappa_grid;
    vector<vector<double> > PDF_grid;
    
    vector<FlatInhomogeneousUniverseLCDM*> universes;
    vector<GalaxySample3D*> galaxy_samples_3D;
    vector<ProjectedGalaxySample*> projected_galaxy_samples;
    
};

Global_Universes global_universes;




cosmological_model set_cosmological_model(double Omega_m, double Omega_b, double Omega_r, double Omega_L, double sigma_8, double n_s, double h_100, double w0, double w1){
  
  cosmological_model cosmo;
  
  cosmo.Omega_m = Omega_m;
  cosmo.Omega_r = Omega_r; //Planck Omega_r
  cosmo.Omega_b = Omega_b;
  cosmo.Omega_L = Omega_L;
  cosmo.Omega_k = 1.0 - cosmo.Omega_m - cosmo.Omega_r - cosmo.Omega_L;
  
  cosmo.n_s = n_s;
  cosmo.h_100 = h_100;
  cosmo.theta_27 = 1.0094444444444444; //Planck theta_27
  cosmo.sigma_8 = sigma_8;
  cosmo.w0 = w0;
  cosmo.w1 = w1;
  
  return cosmo;
    
}


extern "C" void initialise_new_Universe(double a_initial, double a_final, double Omega_m, double Omega_b, double Omega_r, double Omega_L, double sigma_8, double n_s, double h_100, double w0, double w1){
  
  cosmological_model cosmo = set_cosmological_model(Omega_m, Omega_b, Omega_r, Omega_L, sigma_8, n_s, h_100, w0, w1);
  global_universes.cosmo = cosmo;
    
  global_universes.universes.push_back(new FlatInhomogeneousUniverseLCDM(cosmo, a_initial, a_final));
  int N_Universe = global_universes.universes.size();
    
}


extern "C" void add_3D_galaxy_sample(int index_of_universe, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1, double a2, int bias_model){
  BIAS_MODEL b_model = static_cast<BIAS_MODEL>(bias_model);
  global_universes.galaxy_samples_3D.push_back(new GalaxySample3D(global_universes.universes[index_of_universe], z, density_in_Mpc_over_h_cubed, b1, b2, a0, a1, a2, b_model));
}


extern "C" void add_3D_galaxy_sample_from_table(int index_of_universe, double z, double density_in_Mpc_over_h_cubed, const char *bias_table){
  global_universes.galaxy_samples_3D.push_back(new GalaxySample3D(global_universes.universes[index_of_universe], z, density_in_Mpc_over_h_cubed, string(bias_table)));
}

extern "C" int return_N_max_3D(int index_of_galaxy_sample, double R_in_Mpc_over_h, double var_NL_rescale){
  return global_universes.galaxy_samples_3D[index_of_galaxy_sample]->return_N_max_in_3D_tophat(R_in_Mpc_over_h, var_NL_rescale);
}


extern "C" void change_parameters_of_3D_galaxy_sample(int index_of_galaxy_sample, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1, double a2, int bias_model){
  BIAS_MODEL b_model = static_cast<BIAS_MODEL>(bias_model);
  global_universes.galaxy_samples_3D[index_of_galaxy_sample]->set_parameters_3D(z, density_in_Mpc_over_h_cubed, b1, b2, a0, a1, a2, b_model);
}

extern "C" double change_b2_to_minimise_negative_densities_3D(int index_of_galaxy_sample, double R_in_Mpc_over_h, double var_NL_rescale){
  return global_universes.galaxy_samples_3D[index_of_galaxy_sample]->set_b2_to_minimise_negative_densities_in_3D_tophat(R_in_Mpc_over_h, var_NL_rescale);
}

extern "C" void return_CiC_PDF_3D(double* P_of_N, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  vector<double> P_of_N_vector = global_universes.galaxy_samples_3D[index_of_galaxy_sample]->return_CiC_PDF_in_3D_tophat(R_in_Mpc_over_h, f_NL, var_NL_rescale);
  
  int number_of_Ns = P_of_N_vector.size();
  
  for(int i = 0; i < number_of_Ns; i++){
    P_of_N[i] = P_of_N_vector[i];
  }
  
}


extern "C" void add_projected_galaxy_sample(int index_of_universe, const char *n_of_z_file, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, double a2, int bias_model){
  BIAS_MODEL b_model = static_cast<BIAS_MODEL>(bias_model);
  global_universes.projected_galaxy_samples.push_back(new ProjectedGalaxySample(global_universes.universes[index_of_universe], density_in_arcmin_squared, b1, b2, a0, a1, a2, b_model, string(n_of_z_file)));
}

extern "C" void add_projected_galaxy_sample_with_IA(int index_of_universe, const char *n_of_z_file, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, double a2, double A_IA, double alpha_IA, int bias_model){
  BIAS_MODEL b_model = static_cast<BIAS_MODEL>(bias_model);
  global_universes.projected_galaxy_samples.push_back(new ProjectedGalaxySample(global_universes.universes[index_of_universe], density_in_arcmin_squared, b1, b2, a0, a1, a2, A_IA, alpha_IA, b_model, string(n_of_z_file)));
}

extern "C" void add_projected_galaxy_sample_with_IA_and_Alens(int index_of_universe, const char *n_of_z_file, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, double a2, double A_IA, double alpha_IA, double a_lens, int bias_model){
  BIAS_MODEL b_model = static_cast<BIAS_MODEL>(bias_model);
  global_universes.projected_galaxy_samples.push_back(new ProjectedGalaxySample(global_universes.universes[index_of_universe], density_in_arcmin_squared, b1, b2, a0, a1, a2, A_IA, alpha_IA, a_lens, b_model, string(n_of_z_file)));
}


extern "C" void add_projected_galaxy_sample_from_table(int index_of_universe, const char *n_of_z_file, double density_in_arcmin_squared, const char *bias_table){
  global_universes.projected_galaxy_samples.push_back(new ProjectedGalaxySample(global_universes.universes[index_of_universe], density_in_arcmin_squared, string(n_of_z_file), string(bias_table)));
}

extern "C" void add_projected_galaxy_sample_with_IA_from_table(int index_of_universe, const char *n_of_z_file, double density_in_arcmin_squared, double A_IA, double alpha_IA, const char *bias_table){
  global_universes.projected_galaxy_samples.push_back(new ProjectedGalaxySample(global_universes.universes[index_of_universe], density_in_arcmin_squared, A_IA, alpha_IA, string(n_of_z_file), string(bias_table)));
}

extern "C" void add_projected_galaxy_sample_with_IA_and_Alens_from_table(int index_of_universe, const char *n_of_z_file, double density_in_arcmin_squared, double A_IA, double alpha_IA, double a_lens, const char *bias_table){
  global_universes.projected_galaxy_samples.push_back(new ProjectedGalaxySample(global_universes.universes[index_of_universe], density_in_arcmin_squared, A_IA, alpha_IA, a_lens, string(n_of_z_file), string(bias_table)));
}



extern "C" int return_N_max_projected(int index_of_galaxy_sample, double theta, double var_NL_rescale){
  return global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_N_max_in_angular_tophat(theta, var_NL_rescale);
}

extern "C" void change_parameters_of_projected_galaxy_sample(int index_of_galaxy_sample, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, double a2, int bias_model, const char *n_of_z_file){
  BIAS_MODEL b_model = static_cast<BIAS_MODEL>(bias_model);
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->set_parameters_projected(density_in_arcmin_squared, b1, b2, a0, a1, a2, b_model, string(n_of_z_file));
}

extern "C" double change_b2_to_minimise_negative_densities_projected(int index_of_galaxy_sample, double theta, double var_NL_rescale){
  return global_universes.projected_galaxy_samples[index_of_galaxy_sample]->set_b2_to_minimise_negative_densities_in_angular_tophat(theta, var_NL_rescale);
}

extern "C" void return_CiC_PDF_projected(double* P_of_N, double theta, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  cout << "CHECK 1\n";
  vector<double> P_of_N_vector = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_CiC_PDF_in_angular_tophat(theta, f_NL, var_NL_rescale);
  
  cout << "CHECK 2\n";
  int number_of_Ns = P_of_N_vector.size();
  cout << "number_of_Ns = " << number_of_Ns << "\n";
  
  for(int i = 0; i < number_of_Ns; i++){
    P_of_N[i] = P_of_N_vector[i];
  }
  
}

extern "C" void return_CiC_PDF_projected_Gauss(double* P_of_N, double theta, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  vector<double> P_of_N_vector = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_CiC_PDF_in_angular_tophat_Gauss(theta, f_NL, var_NL_rescale);
  
  int number_of_Ns = P_of_N_vector.size();
  
  for(int i = 0; i < number_of_Ns; i++){
    P_of_N[i] = P_of_N_vector[i];
  }
  
}

extern "C" void return_CiC_saddle_point_PDF_projected(double* P_of_N, double theta, double f_NL, double var_NL_rescale, int index_of_galaxy_sample, int next_to_leading){
  
  vector<double> P_of_N_vector = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_CiC_saddle_point_PDF_in_angular_tophat(theta, f_NL, var_NL_rescale, next_to_leading);
  
  int number_of_Ns = P_of_N_vector.size();
  
  for(int i = 0; i < number_of_Ns; i++){
    P_of_N[i] = P_of_N_vector[i];
  }
  
}


/*
 * return_R_in_Mpc_over_h_from_angular_scale
 * 
 * This function calculates what transverse scale R the angular scale theta corresponds to at
 * the mean co-moving distance of galaxy sample index_of_galaxy_sample.
 * 
 */
extern "C" double return_R_in_Mpc_over_h_from_angular_scale(int index_of_galaxy_sample, double theta_in_arcmin){
  vector<double> w_vals;
  vector<double> n_of_w_vals;
  vector<double> lensing_kernel_vals;
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_LOS_data(&w_vals, &n_of_w_vals, &lensing_kernel_vals);
  
  double dw, w, w_mean = 0.0;
  for(int i = 0; i < w_vals.size()-1; i++){
    dw = w_vals[i+1]-w_vals[i];
    w = 0.5*(w_vals[i+1]+w_vals[i]);
    w_mean += w*n_of_w_vals[i]*dw;
  }
  
  return (theta_in_arcmin*constants::arcmin)*(w_mean*constants::c_over_e5);
}



extern "C" void set_primordial_skewness(int index_of_universe, int PNG_modus){
  cout << "Initialising primordial skewnesses from exact integrals.\n";
  cout << "NOTE: this can be very slow!";
  global_universes.universes[index_of_universe]->set_sphere_skewnesses(PNG_modus);
}

extern "C" void set_primordial_skewness_from_eps3_powerlaw_approximation(int index_of_universe, int PNG_modus, double R_0_in_Mpc_over_h){
  cout << "Initialising primordial skewnesses from power law approximation\n";
  cout << "eps3(R) = skewness(R)/sigma(R)^3 = A_eps3*(R/R_0)^n_eps3 .\n";
  global_universes.universes[index_of_universe]->set_sphere_skewnesses_from_eps3_powerlaw_approximation(PNG_modus, R_0_in_Mpc_over_h);
}

// skewnesses stored in "skewness_file" should be for f_NL = 1;
extern "C" void set_primordial_skewness_from_file(int index_of_universe, const char *skewness_file){
  cout << "Initialising primordial skewnesses from file\n";
  cout << skewness_file << " .\n";
  cout << "NOTE: skewnesses in this file should be for f_NL = 1";
  global_universes.universes[index_of_universe]->set_sphere_skewnesses_from_file(skewness_file);
}




extern "C" void set_primordial_skewness_2D(int index_of_universe, int PNG_modus){
  cout << "Initialising primordial skewnesses from exact integrals.\n";
  cout << "NOTE: this can be very slow!";
  global_universes.universes[index_of_universe]->set_cylinder_skewnesses(PNG_modus);
}

extern "C" void set_primordial_skewness_from_eps3_powerlaw_approximation_2D(int index_of_universe, int PNG_modus, double R_0_in_Mpc_over_h){
  cout << "Initialising primordial skewnesses from power law approximation\n";
  cout << "eps3(R) = skewness(R)/sigma(R)^3 = A_eps3*(R/R_0)^n_eps3 .\n";
  global_universes.universes[index_of_universe]->set_cylinder_skewnesses_from_eps3_powerlaw_approximation(PNG_modus, R_0_in_Mpc_over_h);
}

extern "C" void clear_universes(){
  
  for(int i = 0; i < global_universes.universes.size();i++) delete global_universes.universes[i];
  for(int i = 0; i < global_universes.galaxy_samples_3D.size();i++) delete global_universes.galaxy_samples_3D[i];
  for(int i = 0; i < global_universes.projected_galaxy_samples.size();i++) delete global_universes.projected_galaxy_samples[i];
  
  global_universes.universes.clear();
  global_universes.galaxy_samples_3D.clear();
  global_universes.projected_galaxy_samples.clear();
  
}

extern "C" void clear_galaxy_samples(){
  
  for(int i = 0; i < global_universes.galaxy_samples_3D.size();i++) delete global_universes.galaxy_samples_3D[i];
  for(int i = 0; i < global_universes.projected_galaxy_samples.size();i++) delete global_universes.projected_galaxy_samples[i];
  global_universes.galaxy_samples_3D.clear();
  global_universes.projected_galaxy_samples.clear();
  
}


extern "C" void return_background_expansion_to_python(double* t_phys, double* eta, double* a, double* H_conf, double* H_conf_prime, int index_of_universe, int conformal_time_steps){
  
  vector<vector<double> > background_expansion = global_universes.universes[index_of_universe]->return_background_expansion(conformal_time_steps);
  
  int n_column = background_expansion.size();
  
  for(int i = 0; i < conformal_time_steps; i++){
    t_phys[i] = background_expansion[0][i];
    eta[i] = background_expansion[1][i];
    a[i] = background_expansion[2][i];
    H_conf[i] = background_expansion[3][i];
    H_conf_prime[i] = background_expansion[4][i];
  }
  
}


extern "C" void return_linear_growth_to_python(double* eta, double* a, double* D, double* D_prime, int index_of_universe, int conformal_time_steps){
  
  vector<vector<double> > linear_growth_history = global_universes.universes[index_of_universe]->return_linear_growth_history(conformal_time_steps);
  
  int n_column = linear_growth_history.size();
  
  for(int i = 0; i < conformal_time_steps; i++){
    eta[i] = linear_growth_history[0][i];
    a[i]   = linear_growth_history[1][i];
    D[i]   = linear_growth_history[2][i];
    D_prime[i] = linear_growth_history[3][i];
  }
  
}


extern "C" void return_flat_matter_dominated_background_expansion_to_python(double* t_phys, double* eta, double* a, double* H_conf, double* H_conf_prime, int conformal_time_steps){
  
  for(int i = 0; i < conformal_time_steps; i++) FlatHomogeneousUniverseLCDM::expansion_in_flat_matter_dominated_universe(a[i], 1.0, &t_phys[i], &eta[i], &H_conf[i], &H_conf_prime[i]);  
  
}

extern "C" void return_flat_radiation_dominated_background_expansion_to_python(double* t_phys, double* eta, double* a, double* H_conf, double* H_conf_prime, int conformal_time_steps){
  
  for(int i = 0; i < conformal_time_steps; i++) FlatHomogeneousUniverseLCDM::expansion_in_flat_radiation_dominated_universe(a[i], &t_phys[i], &eta[i], &H_conf[i], &H_conf_prime[i]);  
  
}

extern "C" void return_flat_Lambda_dominated_background_expansion_to_python(double* t_phys, double* eta, double* a, double* H_conf, double* H_conf_prime, int conformal_time_steps){
  
  for(int i = 0; i < conformal_time_steps; i++) FlatHomogeneousUniverseLCDM::expansion_in_flat_Lambda_dominated_universe(a[i], &t_phys[i], &eta[i], &H_conf[i], &H_conf_prime[i]);  
  
}

extern "C" int return_n_k(){
  return constants::number_of_k;
}

extern "C" void return_power_spectra(double* k_values, double* P_L, double* P_halofit, int index_of_universe, double z){
  
  double eta = global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z));
  
  vector<vector<double> > power_spectra = global_universes.universes[index_of_universe]->return_power_spectra(eta);
  
  for(int i = 0; i < constants::number_of_k; i++){
    k_values[i] = power_spectra[0][i];
    P_L[i]   = power_spectra[1][i];
    P_halofit[i]   = power_spectra[2][i];
  }
  
}

extern "C" void return_CGF(double* delta_L, double* lambda, double* phi, double* phi_prime, double* phi_prime_prime, double* phi_prime_prime_prime, int* N_lambda_uncollapsed, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > phi_data = global_universes.universes[index_of_universe]->compute_phi_of_lambda_3D(z, R_in_comoving_Mpc/constants::c_over_e5, f_NL, var_NL_rescale);
  
  (*N_lambda_uncollapsed) = phi_data[0].size();
  cout << "N_lambda_uncollapsed = " << (*N_lambda_uncollapsed) << '\n';
  for(int i = 0; i < (*N_lambda_uncollapsed); i++){
    lambda[i] = phi_data[0][i];
    phi[i] = phi_data[1][i];
    phi_prime[i] = phi_data[2][i];
    phi_prime_prime[i] = phi_data[3][i];
    phi_prime_prime_prime[i] = phi_data[4][i];
    delta_L[i] = phi_data[5][i]/(1.0+phi_prime[i]);
  }
  
}

extern "C" void return_CGF_cylinder(double* delta_L, double* lambda, double* phi, double* phi_prime, double* phi_prime_prime, double* phi_prime_prime_prime, int* N_lambda_uncollapsed, double z, double R_in_comoving_Mpc, double L_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  double L = L_in_comoving_Mpc/constants::c_over_e5;
  double R = R_in_comoving_Mpc/constants::c_over_e5;
  double e = global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z));
  
  global_universes.universes[index_of_universe]->set_current_P_NL(e);
  double var_NL_Limber = global_universes.universes[index_of_universe]->variance_of_matter_within_R_NL_2D(R)/L;
  double var_NL_exact = global_universes.universes[index_of_universe]->variance_of_matter_within_R_NL_2D(R,L);
  
  vector<vector<double> > phi_data = global_universes.universes[index_of_universe]->compute_phi_tilde_of_lambda_2D(e, R, f_NL, var_NL_rescale/L*var_NL_exact/var_NL_Limber);
  
  (*N_lambda_uncollapsed) = phi_data[0].size();
  for(int i = 0; i < (*N_lambda_uncollapsed); i++){
    lambda[i] = phi_data[0][i];
    phi[i] = phi_data[1][i];
    phi_prime[i] = phi_data[2][i];
    phi_prime_prime[i] = phi_data[3][i];
    phi_prime_prime_prime[i] = phi_data[4][i];
    delta_L[i] = phi_data[5][i]/(1.0+phi_prime[i]);
  }
  
}

extern "C" void return_reduced_phi_tilde(double* delta_L, double* lambda, double* phi, double* phi_prime, double* phi_prime_prime, double* phi_prime_prime_prime, int* N_lambda_uncollapsed, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  double R = R_in_comoving_Mpc/constants::c_over_e5;
  double e = global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z));
  
  
  global_universes.universes[index_of_universe]->set_current_P_NL(e);
  double var_NL_Limber = global_universes.universes[index_of_universe]->variance_of_matter_within_R_NL_2D(R);
  
  vector<vector<double> > phi_data = global_universes.universes[index_of_universe]->compute_phi_tilde_of_lambda_2D(e, R, f_NL, 1.0/var_NL_Limber);
  
  (*N_lambda_uncollapsed) = phi_data[0].size();
  for(int i = 0; i < (*N_lambda_uncollapsed); i++){
    lambda[i] = phi_data[0][i];
    phi[i] = phi_data[1][i];
    phi_prime[i] = phi_data[2][i];
    phi_prime_prime[i] = phi_data[3][i];
    phi_prime_prime_prime[i] = phi_data[4][i];
    delta_L[i] = phi_data[5][i]/(1.0+phi_prime[i]);
  }
  
}

extern "C" void return_CGF_2D(double* lambda, double* phi, double* phi_prime, double* phi_prime_prime, double* phi_prime_prime_prime, int* N_lambda_uncollapsed, double theta_in_arcmin, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->pointer_to_universe();
  vector<vector<double> > phi_data = pointer_to_universe->return_LOS_integrated_phi_of_lambda(theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale, w_values, kernel_values);
  
  (*N_lambda_uncollapsed) = phi_data[0].size();
  for(int i = 0; i < (*N_lambda_uncollapsed); i++){
    lambda[i] = phi_data[0][i];
    phi[i] = phi_data[1][i];
    phi_prime[i] = phi_data[2][i];
    phi_prime_prime[i] = phi_data[3][i];
    phi_prime_prime_prime[i] = phi_data[4][i];
  }
  
}

extern "C" void return_CGF_lensing(double* lambda, double* phi, double* phi_prime, double* phi_prime_prime, double* phi_prime_prime_prime, int* N_lambda_uncollapsed, double theta_in_arcmin, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->pointer_to_universe();
  vector<vector<double> > phi_data = pointer_to_universe->return_LOS_integrated_phi_of_lambda(theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale, w_values, lensing_kernel_values);
  
  (*N_lambda_uncollapsed) = phi_data[0].size();
  for(int i = 0; i < (*N_lambda_uncollapsed); i++){
    lambda[i] = phi_data[0][i];
    phi[i] = phi_data[1][i];
    phi_prime[i] = phi_data[2][i];
    phi_prime_prime[i] = phi_data[3][i];
    phi_prime_prime_prime[i] = phi_data[4][i];
  }
  
}

extern "C" double return_variance_NL_2D(double theta_in_arcmin, double var_NL_rescale, int index_of_galaxy_sample){
  
  return global_universes.projected_galaxy_samples[index_of_galaxy_sample]->compute_variance_in_angular_tophat(theta_in_arcmin, var_NL_rescale);
  
}

extern "C" void return_PDF(double* delta_values, double* PDF, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > PDF_data = global_universes.universes[index_of_universe]->compute_PDF_3D(z, R_in_comoving_Mpc, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
  }
  
}

extern "C" void return_PDF_Francis(double* delta_values, double* PDF, double z, double R_in_comoving_Mpc, double f_NL, double b1, double b2, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > PDF_data = global_universes.universes[index_of_universe]->compute_PDF_3D_Francis(z, R_in_comoving_Mpc, f_NL, b1, b2, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
  }
  
}

extern "C" double return_variance_L_3D(double z, double R_in_comoving_Mpc, int index_of_universe){
  
  return global_universes.universes[index_of_universe]->return_linear_variance(z, R_in_comoving_Mpc);
  
}

extern "C" double return_variance_NL_3D(double z, double R_in_comoving_Mpc, double var_NL_rescale, int index_of_universe){
  
  return global_universes.universes[index_of_universe]->return_non_linear_variance(z, R_in_comoving_Mpc, var_NL_rescale);
  
}

extern "C" double return_skewness_NL_3D(double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  return global_universes.universes[index_of_universe]->return_3D_skewness(z, R_in_comoving_Mpc, f_NL, var_NL_rescale);
  
}


extern "C" void return_PDF_3D_Lagrangian_bias_terms(double* delta_values, double* PDF, double* bias_term_1, double* bias_term_2, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > PDF_data = global_universes.universes[index_of_universe]->compute_PDF_3D(z, R_in_comoving_Mpc, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
    bias_term_1[d] = PDF_data[2][d];
    bias_term_2[d] = PDF_data[3][d];
    // <delta_g | delta_m > = delta_m + b_1_Lagrange*bias_term_1 + b_2_Lagrange*bias_term_2
  }
  
}

extern "C" void return_PDF_2D_Lagrangian_bias_terms(double* delta_values, double* PDF, double* bias_term_1, double* bias_term_2, double z, double theta_in_arcmin, double L_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  double theta = theta_in_arcmin*constants::pi/180.0/60.0;
  double w = global_universes.universes[index_of_universe]->eta_at_a(1.0) - global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z));
  double R_in_comoving_Mpc = theta*w*constants::c_over_e5;
  cout << "R_in_comoving_Mpc = " << R_in_comoving_Mpc << '\n';
  vector<vector<double> > PDF_data = global_universes.universes[index_of_universe]->compute_PDF_2D(z, R_in_comoving_Mpc, L_in_comoving_Mpc, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
    bias_term_1[d] = PDF_data[2][d];
    bias_term_2[d] = PDF_data[3][d];
    // <delta_g | delta_m > = delta_m + b_1_Lagrange*bias_term_1 + b_2_Lagrange*bias_term_2
  }
  
}

extern "C" void return_PDF_2D_Lagrangian_bias_terms_R(double* delta_values, double* PDF, double* bias_term_1, double* bias_term_2, double* bias_term_3, double z, double R_in_comoving_Mpc, double L_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > PDF_data = global_universes.universes[index_of_universe]->compute_PDF_2D(z, R_in_comoving_Mpc, L_in_comoving_Mpc, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
    bias_term_1[d] = PDF_data[2][d];
    bias_term_2[d] = PDF_data[3][d];
    bias_term_3[d] = PDF_data[4][d];
    // <delta_g | delta_m > = bias_term_1 + b_1_Lagrange*bias_term_2 + b_2_Lagrange*bias_term_3
  }
  
}

extern "C" void return_PDF_2D_Lagrangian_bias_and_variance_terms_R(double* delta_values, double* PDF, double* bias_term_1, double* bias_term_2, double* bias_term_3, double* variance_term_1, double* variance_term_2, double* variance_term_3, double z, double R_in_comoving_Mpc, double L_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > PDF_data = global_universes.universes[index_of_universe]->compute_PDF_2D(z, R_in_comoving_Mpc, L_in_comoving_Mpc, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
    // <delta_g | delta_m > = bias_term_1 + b_1_Lagrange*bias_term_2 + b_2_Lagrange*bias_term_3
    bias_term_1[d] = PDF_data[2][d];
    bias_term_2[d] = PDF_data[3][d];
    bias_term_3[d] = PDF_data[4][d];
    // <delta_g^2 | delta_m > = variance_term_1 + b_1_Lagrange*variance_term_2 + b_1_Lagrange^2*variance_term_3
    variance_term_1[d] = PDF_data[5][d];
    variance_term_2[d] = PDF_data[6][d];
    variance_term_3[d] = PDF_data[7][d];
  }
  
}

extern "C" void return_PDF_at_chosen_deltas_2D_Lagrangian_bias_terms(int N_delta, double* delta_values, double* PDF, double* bias_term_1, double* bias_term_2, double* bias_term_3, double z, double theta_in_arcmin, double L_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  double theta = theta_in_arcmin*constants::pi/180.0/60.0;
  double w = global_universes.universes[index_of_universe]->eta_at_a(1.0) - global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z));
  double R_in_comoving_Mpc = theta*w*constants::c_over_e5;
  cout << "R_in_comoving_Mpc = " << R_in_comoving_Mpc << '\n';
  
  vector<double> delta_vector(N_delta, 0.0);
  for(int d = 0; d < N_delta; d++){
    delta_vector[d] = delta_values[d];
  }
  
  vector<vector<double> > PDF_data = global_universes.universes[index_of_universe]->compute_PDF_at_chosen_deltas_2D(&delta_vector, z, R_in_comoving_Mpc, L_in_comoving_Mpc, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
    bias_term_1[d] = PDF_data[2][d];
    bias_term_2[d] = PDF_data[4][d];
    bias_term_3[d] = PDF_data[5][d];
    // <delta_g | delta_m > = bias_term_1 + b_1_Lagrange*bias_term_2 + b_2_Lagrange*bias_term_3
  }
  
}

extern "C" int return_Ndelta(){
  return constants::N_delta_values_for_PDFs;
}

extern "C" int return_Nlambda(int index_of_universe){
  return global_universes.universes[index_of_universe]->return_N_of_lambda();
}

extern "C" int return_Nlambda_2D(int index_of_universe){
  return global_universes.universes[index_of_universe]->return_N_of_lambda_2D();
}

extern "C" double return_var_NL(double z, double R_in_Mpc_over_h, int index_of_universe){
  // ISSUE: var_NL_rescale non implemented here
  return global_universes.universes[index_of_universe]->return_non_linear_variance(z, R_in_Mpc_over_h, 1.0);
}


extern "C" double return_var_L(double z, double R_in_Mpc_over_h, int index_of_universe){
  return global_universes.universes[index_of_universe]->return_linear_variance(z, R_in_Mpc_over_h);
}

extern "C" void compute_projected_PDF_incl_CMB_kappa(int *Nd, int *Nk, double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, int index_of_galaxy_sample){
  vector<vector<double> > *d_grid = &global_universes.delta_grid;
  vector<vector<double> > *k_grid = &global_universes.kappa_grid;
  vector<vector<double> > *p_grid = &global_universes.PDF_grid;
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->pointer_to_universe()->compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point(theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale, kappa_min, kappa_max, 0.0, w_values, kernel_values, d_grid, k_grid, p_grid);
  (*Nd) = p_grid->size();
  (*Nk) = (*p_grid)[0].size();
}

extern "C" void return_projected_PDF_incl_CMB_kappa(double *delta_grid, double *kappa_grid, double *PDF_grid){
  int Nd = global_universes.PDF_grid.size();
  int Nk = global_universes.PDF_grid[0].size();
  int index = 0;
  for(int k = 0; k < Nk; k++){
    for(int d = 0; d < Nd; d++){
      delta_grid[index] = global_universes.delta_grid[d][k];
      kappa_grid[index] = global_universes.kappa_grid[d][k];
      PDF_grid[index] = global_universes.PDF_grid[d][k];
      index++;
    }
  }
  
}



extern "C" void print_joint_PDF_Ng_kappaCMB(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, int index_of_galaxy_sample){
  
  vector<vector<double> > n_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_joint_saddle_point_PDF_Ng_kappaCMB_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, kappa_min, kappa_max, &n_grid, &k_grid, &p_grid);
  int NN = n_grid.size();
  int Nk = n_grid[0].size();
  int index = 0;
  
  cout << NN << '\n';
  cout << Nk << '\n';
  
  FILE *F = fopen("PDF_test", "w");
  fclose(F);
  F = fopen("N_test", "w");
  fclose(F);
  F = fopen("kappa_test", "w");
  fclose(F);
  fstream out_n, out_k, out_p;
  out_n.open("N_test");
  out_k.open("kappa_test");
  out_p.open("PDF_test");
  out_n << scientific << setprecision(5);
  out_k << scientific << setprecision(5);
  out_p << scientific << setprecision(5);
  
  for(int n = 0; n < NN; n++){
    for(int k = 0; k < Nk; k++){
      out_n << n_grid[n][k] << "  ";
      out_k << k_grid[n][k] << "  ";
      out_p << p_grid[n][k] << "  ";
    }
    out_n << '\n';
    out_k << '\n';
    out_p << '\n';
  }
  
  out_n.close();
  out_k.close();
  out_p.close();
  
}

extern "C" void print_joint_PDF_Ng_kappaCMB_noisy(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_CMB_noise_variance, int index_of_galaxy_sample){
  
  vector<vector<double> > n_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_joint_saddle_point_PDF_Ng_kappaCMB_noisy_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, kappa_min, kappa_max, kappa_CMB_noise_variance, &n_grid, &k_grid, &p_grid);
  int NN = n_grid.size();
  int Nk = n_grid[0].size();
  int index = 0;
  
  cout << NN << '\n';
  cout << Nk << '\n';
  
  FILE *F = fopen("PDF_noisyKappa_test", "w");
  fclose(F);
  F = fopen("N_noisyKappa_test", "w");
  fclose(F);
  F = fopen("kappa_noisyKappa_test", "w");
  fclose(F);
  fstream out_n, out_k, out_p;
  out_n.open("N_noisyKappa_test");
  out_k.open("kappa_noisyKappa_test");
  out_p.open("PDF_noisyKappa_test");
  out_n << scientific << setprecision(5);
  out_k << scientific << setprecision(5);
  out_p << scientific << setprecision(5);
  
  for(int n = 0; n < NN; n++){
    for(int k = 0; k < Nk; k++){
      out_n << n_grid[n][k] << "  ";
      out_k << k_grid[n][k] << "  ";
      out_p << p_grid[n][k] << "  ";
    }
    out_n << '\n';
    out_k << '\n';
    out_p << '\n';
  }
  
  out_n.close();
  out_k.close();
  out_p.close();
  
}


extern "C" void return_joint_PDF_Ng_kappaCMB_noisy(double* joint_rebinned_PDF, double theta_in_arcmin, double f_NL, double var_NL_rescale, double N_min, double N_max, double kappa_min, double kappa_max, double kappa_CMB_noise_variance, int N_bin, int index_of_galaxy_sample){
  
  vector<vector<double> > n_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  
  // Need to use global kappa range here since convolution with kappa noise has to be done with the full PDF
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_joint_saddle_point_PDF_Ng_kappaCMB_noisy_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, constants::kappa_CMB_min, constants::kappa_CMB_max, kappa_CMB_noise_variance, &n_grid, &k_grid, &p_grid);
  
  int NN = n_grid.size();
  int Nk = n_grid[0].size();
  
  //double dk = (constants::kappa_CMB_max - constants::kappa_CMB_min)/double(Nk-1);
  double dk = k_grid[0][1] - k_grid[0][0];
  double dN = 1.0;
  double dV = dk*dN;
  
  double dk_rebin = (kappa_max - kappa_min)/double(N_bin);
  double dN_rebin = (N_max - N_min)/double(N_bin);
  double dV_rebin = dk_rebin*dN_rebin;
  
  vector<vector<double> > PDF_grid_rebin_N(N_bin, vector<double>(Nk,0.0));
  int n_rebin;
  double N;
  double N_bin_min, N_bin_max;
  double norm = 0.0;
  double norm2 = 0.0;
  for(int k = 0; k < Nk; k++){
    for(int n = 0; n < NN; n++){
      N = n_grid[n][k];
      n_rebin = int((N - N_min)/dN_rebin);
      if(n_rebin > -1 && n_rebin < N_bin && N >= N_min){
        PDF_grid_rebin_N[n_rebin][k] += p_grid[n][k];
        norm += p_grid[n][k]*dk*dN;
      }
      norm2 += p_grid[n][k]*dk*dN;
    }
    
  }
  
  cout << "Testing the norm again: " << norm << '\n';
  cout << "Testing the norm again: " << norm2 << '\n';
  cout << "N_min: " << N_min << '\n';
  cout << "N_bin: " << N_bin << '\n';
  cout << "dN_rebin: " << dN_rebin << '\n';
  
  vector<vector<double> > PDF_grid_rebinned(N_bin, vector<double>(N_bin,0.0));
  int k_rebin;
  double kappa;
  double kappa_bin_min, kappa_bin_max;
  
  norm = 0.0;
  for(int n = 0; n < N_bin; n++){
    
    for(int k = 0; k < Nk; k++){
      kappa = k_grid[0][k];
      k_rebin = int((kappa - kappa_min)/dk_rebin);
      if(k_rebin > -1 && k_rebin < N_bin){
        kappa_bin_min = kappa_min + double(k_rebin)*dk_rebin;
        kappa_bin_max = kappa_min + double(k_rebin+1)*dk_rebin;
        if(kappa >= kappa_bin_min + 0.5*dk && kappa < kappa_bin_max - 0.5*dk){
          // if finer bin is fully contained in coarser bin:
          PDF_grid_rebinned[n][k_rebin] += dk*PDF_grid_rebin_N[n][k];
        }
        else if(kappa >= kappa_bin_min && kappa < kappa_bin_min + 0.5*dk){
          // if finer bin overlaps with coarser bin from below:
          PDF_grid_rebinned[n][k_rebin] += (kappa+0.5*dk-kappa_bin_min)*PDF_grid_rebin_N[n][k];
          if(k_rebin>0) PDF_grid_rebinned[n][k_rebin-1] += (kappa_bin_min - kappa + 0.5*dk)*PDF_grid_rebin_N[n][k];
        }
        else if(kappa >= kappa_bin_max - 0.5*dk && kappa < kappa_bin_max ){
          // if finer bin overlaps with coarser bin from above:
          PDF_grid_rebinned[n][k_rebin] += (kappa_bin_max-kappa+0.5*dk)*PDF_grid_rebin_N[n][k];
          if(k_rebin<N_bin-1) PDF_grid_rebinned[n][k_rebin+1] += (kappa + 0.5*dk - kappa_bin_max)*PDF_grid_rebin_N[n][k];
        }
      }
      else if(k_rebin == -1 && kappa > kappa_min - 0.5*dk){
        PDF_grid_rebinned[n][0] += (kappa + 0.5*dk - kappa_min)*PDF_grid_rebin_N[n][k];
      }
      else if(k_rebin == N_bin && kappa < kappa_max + 0.5*dk){
        PDF_grid_rebinned[n][N_bin-1] += (kappa_max - kappa + 0.5*dk)*PDF_grid_rebin_N[n][k];
      }
      
      norm += PDF_grid_rebin_N[n][k]*dk;
    }
    
    for(int k = 0; k < N_bin; k++){
      PDF_grid_rebinned[n][k] /= dV_rebin;
    }
    
  }
  
  int index = 0;
  norm = 0.0;
  for(int n = 0; n < N_bin; n++){
    for(int k = 0; k < N_bin; k++){
      joint_rebinned_PDF[index] = PDF_grid_rebinned[n][k];
      norm += joint_rebinned_PDF[index];
      index++;
    }
  }
  
  norm *= dV_rebin;
  cout << "Norm is now: " << norm << '\n';
  
}


/*
 * CURRENTLY WORKING HERE
 * -->
 * 
 */


extern "C" void return_joint_PDF_Ng_kappa_noisy(double* joint_PDF, double* N_values, double* kappa_values, double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_noise_variance, int index_of_lens_sample, int index_of_source_sample){
  
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  
  global_universes.projected_galaxy_samples[index_of_source_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  
  
  vector<vector<double> > n_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  
  global_universes.projected_galaxy_samples[index_of_lens_sample]->return_joint_saddle_point_PDF_Ng_kappa_noisy_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, -1, +1, kappa_noise_variance, w_values, lensing_kernel_values, &n_grid, &k_grid, &p_grid);
  
  int N_N = p_grid.size();
  int N_kappa = p_grid[0].size();
  
  for(int iN = 0; iN < N_N; iN++){
    N_values[iN] = n_grid[iN][0];
  }
  
  for(int ik = 0; ik < N_kappa; ik++){
    kappa_values[ik] = k_grid[0][ik];
  }
  
  int index = 0;
  for(int iN = 0; iN < N_N; iN++){
    for(int ik = 0; ik < N_kappa; ik++){
      joint_PDF[index] = p_grid[iN][ik];
      index++;
    }
  }
  
}


extern "C" void return_joint_rebinned_PDF_Ng_kappa_noisy(double* joint_rebinned_PDF, double theta_in_arcmin, double f_NL, double var_NL_rescale, double N_min, double N_max, double kappa_min, double kappa_max, double kappa_noise_variance, int N_bin, int index_of_lens_sample, int index_of_source_sample){
//extern "C" void print_joint_PDF_Ng_kappa_noisy(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_noise_variance, int index_of_lens_sample, int index_of_source_sample){
  
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  
  global_universes.projected_galaxy_samples[index_of_source_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  
  
  vector<vector<double> > n_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  
  global_universes.projected_galaxy_samples[index_of_lens_sample]->return_joint_saddle_point_PDF_Ng_kappa_noisy_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, kappa_min, kappa_max, kappa_noise_variance, w_values, lensing_kernel_values, &n_grid, &k_grid, &p_grid);
  
  //vector<vector<double> > n_grid;
  //vector<vector<double> > k_grid;
  //vector<vector<double> > p_grid;
  
  // Need to use global kappa range here since convolution with kappa noise has to be done with the full PDF
  //global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_joint_saddle_point_PDF_Ng_kappaCMB_noisy_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, constants::kappa_CMB_min, constants::kappa_CMB_max, kappa_CMB_noise_variance, &n_grid, &k_grid, &p_grid);
  
  int NN = n_grid.size();
  int Nk = n_grid[0].size();
  
  //double dk = (constants::kappa_CMB_max - constants::kappa_CMB_min)/double(Nk-1);
  double dk = (kappa_max - kappa_min)/double(Nk-1);
  double dN = 1.0;
  double dV = dk*dN;
  
  double dk_rebin = (kappa_max - kappa_min)/double(N_bin);
  double dN_rebin = (N_max - N_min)/double(N_bin);
  double dV_rebin = dk_rebin*dN_rebin;
  
  vector<vector<double> > PDF_grid_rebin_N(N_bin, vector<double>(Nk,0.0));
  int n_rebin;
  double N;
  double N_bin_min, N_bin_max;
  double norm = 0.0;
  for(int k = 0; k < Nk; k++){
    for(int n = 0; n < NN; n++){
      N = n_grid[n][k];
      n_rebin = int((N - N_min)/dN_rebin);
      if(n_rebin > -1 && n_rebin < N_bin && N >= N_min){
        PDF_grid_rebin_N[n_rebin][k] += p_grid[n][k];
        norm += p_grid[n][k]*dk*dN;
      }
    }
    
  }
  
  vector<vector<double> > PDF_grid_rebinned(N_bin, vector<double>(N_bin,0.0));
  int k_rebin;
  double kappa;
  double kappa_bin_min, kappa_bin_max;
  
  norm = 0.0;
  for(int n = 0; n < N_bin; n++){
    
    for(int k = 0; k < Nk; k++){
      kappa = k_grid[0][k];
      k_rebin = int((kappa - kappa_min)/dk_rebin);
      if(k_rebin > -1 && k_rebin < N_bin){
        kappa_bin_min = kappa_min + double(k_rebin)*dk_rebin;
        kappa_bin_max = kappa_min + double(k_rebin+1)*dk_rebin;
        if(kappa >= kappa_bin_min + 0.5*dk && kappa < kappa_bin_max - 0.5*dk){
          // if finer bin is fully contained in coarser bin:
          PDF_grid_rebinned[n][k_rebin] += dk*PDF_grid_rebin_N[n][k];
        }
        else if(kappa >= kappa_bin_min && kappa < kappa_bin_min + 0.5*dk){
          // if finer bin overlaps with coarser bin from below:
          PDF_grid_rebinned[n][k_rebin] += (kappa+0.5*dk-kappa_bin_min)*PDF_grid_rebin_N[n][k];
          if(k_rebin>0) PDF_grid_rebinned[n][k_rebin-1] += (kappa_bin_min - kappa + 0.5*dk)*PDF_grid_rebin_N[n][k];
        }
        else if(kappa >= kappa_bin_max - 0.5*dk && kappa < kappa_bin_max ){
          // if finer bin overlaps with coarser bin from above:
          PDF_grid_rebinned[n][k_rebin] += (kappa_bin_max-kappa+0.5*dk)*PDF_grid_rebin_N[n][k];
          if(k_rebin<N_bin-1) PDF_grid_rebinned[n][k_rebin+1] += (kappa + 0.5*dk - kappa_bin_max)*PDF_grid_rebin_N[n][k];
        }
      }
      else if(k_rebin == -1 && kappa > kappa_min - 0.5*dk){
        PDF_grid_rebinned[n][0] += (kappa + 0.5*dk - kappa_min)*PDF_grid_rebin_N[n][k];
      }
      else if(k_rebin == N_bin && kappa < kappa_max + 0.5*dk){
        PDF_grid_rebinned[n][N_bin-1] += (kappa_max - kappa + 0.5*dk)*PDF_grid_rebin_N[n][k];
      }
      
      norm += PDF_grid_rebin_N[n][k]*dk;
    }
    
    for(int k = 0; k < N_bin; k++){
      PDF_grid_rebinned[n][k] /= dV_rebin;
    }
    
  }
  
  int index = 0;
  for(int n = 0; n < N_bin; n++){
    for(int k = 0; k < N_bin; k++){
      joint_rebinned_PDF[index] = PDF_grid_rebinned[n][k];
      index++;
    }
  }
  
}




/*
 * <--
 * CURRENTLY WORKING HERE
 * 
 */


extern "C" void print_joint_PDF_Ng_kappa_noisy(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_noise_variance, int index_of_lens_sample, int index_of_source_sample){
  
  
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  
  global_universes.projected_galaxy_samples[index_of_source_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  
  
  vector<vector<double> > n_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  
  global_universes.projected_galaxy_samples[index_of_lens_sample]->return_joint_saddle_point_PDF_Ng_kappa_noisy_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, kappa_min, kappa_max, kappa_noise_variance, w_values, lensing_kernel_values, &n_grid, &k_grid, &p_grid);
  int NN = n_grid.size();
  int Nk = n_grid[0].size();
  int index = 0;
  
  cout << NN << '\n';
  cout << Nk << '\n';
  
  FILE *F = fopen("PDF_noisyKappa_test", "w");
  fclose(F);
  F = fopen("N_noisyKappa_test", "w");
  fclose(F);
  F = fopen("kappa_noisyKappa_test", "w");
  fclose(F);
  fstream out_n, out_k, out_p;
  out_n.open("N_noisyKappa_test");
  out_k.open("kappa_noisyKappa_test");
  out_p.open("PDF_noisyKappa_test");
  out_n << scientific << setprecision(5);
  out_k << scientific << setprecision(5);
  out_p << scientific << setprecision(5);
  
  for(int n = 0; n < NN; n++){
    for(int k = 0; k < Nk; k++){
      out_n << n_grid[n][k] << "  ";
      out_k << k_grid[n][k] << "  ";
      out_p << p_grid[n][k] << "  ";
    }
    out_n << '\n';
    out_k << '\n';
    out_p << '\n';
  }
  
  out_n.close();
  out_k.close();
  out_p.close();
  
}

/*
extern "C" void return_joint_PDF_Ng_kappaCMB_noisy_Takahashi(double* joint_rebinned_PDF, double theta_in_arcmin, double f_NL, double var_NL_rescale, double N_min, double N_max, double kappa_min, double kappa_max, double kappa_noise_variance, int N_bin, int index_of_lens_sample, int index_of_source_sample){
  
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  
  global_universes.projected_galaxy_samples[index_of_source_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  
  vector<vector<double> > n_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  
  // Need to use global kappa range here since convolution with kappa noise has to be done with the full PDF
  global_universes.projected_galaxy_samples[index_of_lens_sample]->return_joint_saddle_point_PDF_Ng_kappa_noisy_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, constants::kappa_CMB_min, constants::kappa_CMB_max, kappa_CMB_noise_variance, w_values, lensing_kernel_values, &n_grid, &k_grid, &p_grid);
  //global_universes.projected_galaxy_samples[index_of_lens_sample]->return_joint_saddle_point_PDF_Ng_kappa_noisy_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, kappa_min, kappa_max, kappa_noise_variance, w_values, lensing_kernel_values, &n_grid, &k_grid, &p_grid);
  
  int NN = n_grid.size();
  int Nk = n_grid[0].size();
  
  double dk = (constants::kappa_CMB_max - constants::kappa_CMB_min)/double(Nk-1);
  double dN = 1.0;
  double dV = dk*dN;
  
  double dk_rebin = (kappa_max - kappa_min)/double(N_bin);
  double dN_rebin = (N_max - N_min)/double(N_bin);
  double dV_rebin = dk_rebin*dN_rebin;
  
  vector<vector<double> > PDF_grid_rebin_N(N_bin, vector<double>(Nk,0.0));
  int n_rebin;
  double N;
  double N_bin_min, N_bin_max;
  double norm = 0.0;
  for(int k = 0; k < Nk; k++){
    for(int n = 0; n < NN; n++){
      N = n_grid[n][k];
      n_rebin = int((N - N_min)/dN_rebin);
      if(n_rebin > -1 && n_rebin < N_bin && N >= N_min){
        PDF_grid_rebin_N[n_rebin][k] += p_grid[n][k];
        norm += p_grid[n][k]*dk*dN;
      }
    }
    
  }
  
  vector<vector<double> > PDF_grid_rebinned(N_bin, vector<double>(N_bin,0.0));
  int k_rebin;
  double kappa;
  double kappa_bin_min, kappa_bin_max;
  
  norm = 0.0;
  for(int n = 0; n < N_bin; n++){
    
    for(int k = 0; k < Nk; k++){
      kappa = k_grid[0][k];
      k_rebin = int((kappa - kappa_min)/dk_rebin);
      if(k_rebin > -1 && k_rebin < N_bin){
        kappa_bin_min = kappa_min + double(k_rebin)*dk_rebin;
        kappa_bin_max = kappa_min + double(k_rebin+1)*dk_rebin;
        if(kappa >= kappa_bin_min + 0.5*dk && kappa < kappa_bin_max - 0.5*dk){
          // if finer bin is fully contained in coarser bin:
          PDF_grid_rebinned[n][k_rebin] += dk*PDF_grid_rebin_N[n][k];
        }
        else if(kappa >= kappa_bin_min && kappa < kappa_bin_min + 0.5*dk){
          // if finer bin overlaps with coarser bin from below:
          PDF_grid_rebinned[n][k_rebin] += (kappa+0.5*dk-kappa_bin_min)*PDF_grid_rebin_N[n][k];
          if(k_rebin>0) PDF_grid_rebinned[n][k_rebin-1] += (kappa_bin_min - kappa + 0.5*dk)*PDF_grid_rebin_N[n][k];
        }
        else if(kappa >= kappa_bin_max - 0.5*dk && kappa < kappa_bin_max ){
          // if finer bin overlaps with coarser bin from above:
          PDF_grid_rebinned[n][k_rebin] += (kappa_bin_max-kappa+0.5*dk)*PDF_grid_rebin_N[n][k];
          if(k_rebin<N_bin-1) PDF_grid_rebinned[n][k_rebin+1] += (kappa + 0.5*dk - kappa_bin_max)*PDF_grid_rebin_N[n][k];
        }
      }
      else if(k_rebin == -1 && kappa > kappa_min - 0.5*dk){
        PDF_grid_rebinned[n][0] += (kappa + 0.5*dk - kappa_min)*PDF_grid_rebin_N[n][k];
      }
      else if(k_rebin == N_bin && kappa < kappa_max + 0.5*dk){
        PDF_grid_rebinned[n][N_bin-1] += (kappa_max - kappa + 0.5*dk)*PDF_grid_rebin_N[n][k];
      }
      
      norm += PDF_grid_rebin_N[n][k]*dk;
    }
    
    for(int k = 0; k < N_bin; k++){
      PDF_grid_rebinned[n][k] /= dV_rebin;
    }
    
  }
  
  int index = 0;
  for(int n = 0; n < N_bin; n++){
    for(int k = 0; k < N_bin; k++){
      joint_rebinned_PDF[index] = PDF_grid_rebinned[n][k];
      index++;
    }
  }
  
}
*/


extern "C" void return_convergence_PDF_from_single_z(double* kappa_values, double* PDF, double z_source, double theta_in_arcmin, double f_NL, double var_NL_rescale, int index_of_universe){
  
  
  double eta_0 = global_universes.universes[index_of_universe]->eta_at_a(1.0);
  double eta_source = global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z_source));
  double w_source = eta_0-eta_source;
  double w = 0.0;
  double dw, eta_1, eta_2, scale_1, scale_2;
  
  vector<double> w_boundaries(0, 0.0);
  
  while(w < w_source){
    w_boundaries.push_back(w);
    eta_1 = eta_0-w;
    scale_1 = global_universes.universes[index_of_universe]->a_at_eta(eta_1);
    scale_2 = scale_1 - constants::maximal_da;
    eta_2 = global_universes.universes[index_of_universe]->eta_at_a(scale_2);
    w += min(constants::maximal_dw, eta_1 - eta_2);
  }
  w_boundaries.push_back(w_source);
  
  int N_w = w_boundaries.size();
  vector<double> lensing_kernel_values(N_w, 0.0);
  
  for(int i = 0; i < N_w-1; i++){
    w = 0.5*(w_boundaries[i]+w_boundaries[i+1]);
    scale_1 = global_universes.universes[index_of_universe]->a_at_eta(eta_0-w);
    lensing_kernel_values[i] = 1.5/w_source*global_universes.universes[index_of_universe]->return_Omega_m();
    lensing_kernel_values[i] *= w*(w_source-w)/scale_1;
  }
  
  vector<vector<double> > PDF_data = global_universes.universes[index_of_universe]->compute_LOS_projected_PDF(w_boundaries, lensing_kernel_values, theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    kappa_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
  }
  
}

extern "C" void return_convergence_PDF_incl_galaxies(double* kappa_values, double* PDF, double* expec_dgal_given_kappa, double theta_in_arcmin, double f_NL, double var_NL_rescale, double bL_1, double bL_2, int index_of_source_sample, int index_of_lens_sample){
  
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_source_sample]->pointer_to_universe();

  
  vector<double> w_values_sample_1, w_values_sample_2;
  vector<double> kernel_values_sample_1, kernel_values_sample_2;
  vector<double> lensing_kernel_values_sample_1, lensing_kernel_values_sample_2;
  global_universes.projected_galaxy_samples[index_of_source_sample]->return_LOS_data(&w_values_sample_1, &kernel_values_sample_1, &lensing_kernel_values_sample_1);
  global_universes.projected_galaxy_samples[index_of_lens_sample]->return_LOS_data(&w_values_sample_2, &kernel_values_sample_2, &lensing_kernel_values_sample_2);
  
  vector<vector<double> > rebinned_data = return_joint_binning(w_values_sample_1, lensing_kernel_values_sample_1, w_values_sample_2, kernel_values_sample_2);
  
  vector<double> w_values = rebinned_data[0];
  vector<double> kernel_values_sources = rebinned_data[1];
  vector<double> kernel_values_lenses = rebinned_data[2];
  
  vector<vector<double> > PDF_data = pointer_to_universe->compute_LOS_projected_PDF_incl_galaxies(w_values, kernel_values_sources, kernel_values_lenses, theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale);
  for(int d = 0; d < PDF_data[0].size(); d++){
    kappa_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
    expec_dgal_given_kappa[d] = PDF_data[4][d] + bL_1*PDF_data[2][d] + bL_2*PDF_data[3][d];
  }
  
}

extern "C" void return_convergence_PDF_from_source_sample(double* kappa_values, double* PDF, double theta_in_arcmin, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  vector<vector<double> > PDF_data = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_kappa_PDF_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    kappa_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
  }
  
}

extern "C" void return_matter_PDF_from_tracer_sample(double* delta_values, double* PDF, double theta_in_arcmin, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  vector<vector<double> > PDF_data = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_matter_PDF_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
  }
  
}

extern "C" void return_matter_saddle_point_PDF_from_tracer_sample(double* delta_values, double* PDF, double theta_in_arcmin, double f_NL, double var_NL_rescale, int next_to_leading, int index_of_galaxy_sample){
  
  vector<vector<double> > PDF_data = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_matter_saddle_point_PDF_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale, next_to_leading);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
  }
  
}

extern "C" void configure_FLASK_for_delta_g_and_CMB_kappa(int l_max, double theta_in_arcmin, double bias, double r, int index_of_galaxy_sample){
  
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->pointer_to_universe();

  
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  
  int n_time = w_values.size()-1;
  double a, eta_0, w, dw, w_last_scattering;
  eta_0 = pointer_to_universe->eta_at_a(1.0);
  w_last_scattering = eta_0-pointer_to_universe->eta_at_a(1.0/(1.0+constants::z_last_scattering));
  
  // need to extend projection kernel to surface of last scattering:
  vector<double> w_values_extended_to_last_scattering = w_values;
  vector<double> kernel_values_extended_to_last_scattering = kernel_values;
  dw = min(constants::maximal_dw, w_last_scattering - w_values_extended_to_last_scattering[n_time]);
  w = w_values_extended_to_last_scattering[n_time] + dw;
  while(w < w_last_scattering){
    w_values_extended_to_last_scattering.push_back(w);
    kernel_values_extended_to_last_scattering.push_back(0.0);
    n_time = w_values_extended_to_last_scattering.size()-1;
    dw = min(constants::maximal_dw, w_last_scattering - w);
    w = w + dw;
  }
  w_values_extended_to_last_scattering.push_back(w);
  kernel_values_extended_to_last_scattering.push_back(0.0);
  n_time = w_values_extended_to_last_scattering.size()-1;
  
  vector<double> w_values_bin_center(n_time, 0.0);
  for(int i = 0; i < n_time; i++){
    w_values_bin_center[i] = 0.5*(w_values_extended_to_last_scattering[i+1]+w_values_extended_to_last_scattering[i]);
  }
  
  vector<double> CMB_lensing_kernel(n_time, 0.0);
  vector<double> CMB_lensing_kernel_overlap(n_time, 0.0);
  vector<double> CMB_lensing_kernel_nonoverlap(n_time, 0.0);
  double w_min = w_values_extended_to_last_scattering[0];
  double w_max = w_values_extended_to_last_scattering[n_time];
  double norm = 0.0;
  for(int i = 0; i < n_time; i++){
    w = w_values_bin_center[i];
    a = pointer_to_universe->a_at_eta(eta_0-w);
    CMB_lensing_kernel[i] = 1.5*pointer_to_universe->return_Omega_m()*w*(w_last_scattering-w)/w_last_scattering/a;
    if(kernel_values_extended_to_last_scattering[i] > 0.0){
      CMB_lensing_kernel_overlap[i] = CMB_lensing_kernel[i];
    }
    else{
      CMB_lensing_kernel_nonoverlap[i] = CMB_lensing_kernel[i];
    }
      
  }
  
  
  vector<vector<double> > joint_kernel_values(4, vector<double>(n_time, 0.0));
  
  joint_kernel_values[0] = kernel_values_extended_to_last_scattering;
  joint_kernel_values[1] = CMB_lensing_kernel;
  joint_kernel_values[2] = CMB_lensing_kernel_overlap;
  joint_kernel_values[3] = CMB_lensing_kernel_nonoverlap;
  
  
  vector<vector<vector<double> > > C_ells = pointer_to_universe->return_LOS_integrated_C_ells(l_max, w_values_extended_to_last_scattering, joint_kernel_values);
  
  vector<vector<double> > second_moments = pointer_to_universe->return_LOS_integrated_2nd_moments(theta_in_arcmin*constants::arcmin, 0.0, 1.0, w_values_extended_to_last_scattering, joint_kernel_values);
  
  vector<vector<vector<double> > > third_moments = pointer_to_universe->return_LOS_integrated_3rd_moments(theta_in_arcmin*constants::arcmin, 0.0, 1.0, w_values_extended_to_last_scattering, joint_kernel_values);
  
  int N_ell = C_ells[0][0].size();
  
  
  double lambda_m = lognormal_tools::get_delta0(second_moments[0][0], third_moments[0][0][0]);
  double lambda_k_overlap_correlated = lognormal_tools::get_kappa0(lambda_m, second_moments[0][0], second_moments[0][1], third_moments[0][0][1]);
  
  double var_kappa_overlap = second_moments[2][2];
  double var_kappa_overlap_correlated = lognormal_tools::get_var_kappa(lambda_m, lambda_k_overlap_correlated, second_moments[0][1], third_moments[0][1][1]);
  double var_kappa_overlap_uncorrelated = var_kappa_overlap - var_kappa_overlap_correlated;

  double skewness_kappa_overlap = third_moments[2][2][2];
  double skewness_kappa_overlap_correlated = lognormal_tools::get_skew_from_delta_0(lambda_k_overlap_correlated, var_kappa_overlap_correlated);
  double skewness_kappa_overlap_uncorrelated = skewness_kappa_overlap - skewness_kappa_overlap_correlated;
  
  double var_kappa_uncorrelated = second_moments[3][3] + var_kappa_overlap_uncorrelated;
  double skewness_uncorrelated = third_moments[3][3][3] + skewness_kappa_overlap_uncorrelated;
  double lambda_k_uncorrelated = lognormal_tools::get_delta0(var_kappa_uncorrelated, skewness_uncorrelated);
  
  
  FILE *F = fopen("FLASK_config_Cells", "w");
  fclose(F);
  fstream out;
  out.open("FLASK_config_Cells");
  out << scientific << setprecision(5);
  
  out << "# variances: var_g = ";
  out << second_moments[0][0]*bias*bias << " ;  cov_gk = ";
  out << second_moments[0][1]*bias*r << " ;  var_k_correlated = ";
  out << var_kappa_overlap_correlated << " ;  var_k_uncorrelated = ";
  out << var_kappa_uncorrelated << "\n";
  
  out << "# shift params: lambda_g = ";
  out << lambda_m*bias << " ;  lambda_k_correlated = ";
  out << lambda_k_overlap_correlated << " ;  lambda_k_uncorrelated = ";
  out << lambda_k_uncorrelated << "\n";
  
  out << "# ell          C_gg(ell)          C_gk(ell)          C_kk_correlated(ell)          C_kk_uncorrelated(ell)\n";
  
  for(int l = 0; l < N_ell; l++){
    out << l << setw(15);
    out << C_ells[0][0][l]*bias*bias << setw(15);
    out << C_ells[0][1][l]*bias*r << setw(15);
    out << C_ells[2][2][l]*var_kappa_overlap_correlated/var_kappa_overlap << setw(15);
    out << C_ells[3][3][l] + C_ells[2][2][l]*(1.0-var_kappa_overlap_correlated/var_kappa_overlap) << '\n';
  }
  
  out.close();
  
  
}



extern "C" double return_lognormal_shift_for_individual_FLASK_bin(double theta_in_arcmin, int index_of_galaxy_sample, int sources0_or_lenses1){
  
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->pointer_to_universe();
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel;
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel);
  
  int n_time = w_values.size();
  
  for(int i = 0; i < n_time; i++){
    cout << i << "  " << w_values[i] << "  " << kernel_values[i] << "  " << lensing_kernel[i] << '\n';
  }
  
  vector<vector<double> > joint_kernel_values(1, vector<double>(n_time, 0.0));
  if(sources0_or_lenses1 == 0){
    joint_kernel_values[0] = lensing_kernel;
  }
  else if(sources0_or_lenses1 == 1){
    joint_kernel_values[0] = kernel_values;
  }
  else{
    error_handling::general_error_message("Argument sources0_or_lenses1 in function return_lognormal_shift_for_individual_FLASK_bin shoud be either 0 or 1.");
  }
  
  vector<vector<double> > second_moments = pointer_to_universe->return_LOS_integrated_2nd_moments(theta_in_arcmin*constants::arcmin, 0.0, 1.0, w_values, joint_kernel_values);
  vector<vector<vector<double> > > third_moments = pointer_to_universe->return_LOS_integrated_3rd_moments(theta_in_arcmin*constants::arcmin, 0.0, 1.0, w_values, joint_kernel_values);
    
  cout << "variance = " << second_moments[0][0] << '\n';
  cout << "skewness = " << third_moments[0][0][0] << '\n';
  
  /***** FLASK configs for sample 1: *****/
  double lambda = lognormal_tools::get_delta0(second_moments[0][0], third_moments[0][0][0]);
  
  return lambda;
  
}



extern "C" void configure_FLASK_for_delta_g_and_CMB_kappa_2_uncorrelated_samples(int l_max, double theta_in_arcmin, double bias_sample_1, double r_sample_1, double bias_sample_2, double r_sample_2, int index_of_galaxy_sample_1, int index_of_galaxy_sample_2){
  
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_galaxy_sample_1]->pointer_to_universe();

  
  vector<double> w_values_sample_1, w_values_sample_2;
  vector<double> kernel_values_sample_1, kernel_values_sample_2;
  vector<double> lensing_kernel_values_sample_1, lensing_kernel_values_sample_2;
  global_universes.projected_galaxy_samples[index_of_galaxy_sample_1]->return_LOS_data(&w_values_sample_1, &kernel_values_sample_1, &lensing_kernel_values_sample_1);
  global_universes.projected_galaxy_samples[index_of_galaxy_sample_2]->return_LOS_data(&w_values_sample_2, &kernel_values_sample_2, &lensing_kernel_values_sample_2);
  
  vector<vector<double> > rebinned_data = return_joint_binning(w_values_sample_1, kernel_values_sample_1, w_values_sample_2, kernel_values_sample_2);
  
  vector<double> w_values = rebinned_data[0];
  kernel_values_sample_1 = rebinned_data[1];
  kernel_values_sample_2 = rebinned_data[2];
  
  
  int n_time = w_values.size()-1;
  double a, eta_0, w, dw, w_last_scattering;
  eta_0 = pointer_to_universe->eta_at_a(1.0);
  w_last_scattering = eta_0-pointer_to_universe->eta_at_a(1.0/(1.0+constants::z_last_scattering));
  
  // need to extend projection kernel to surface of last scattering:
  vector<double> w_values_extended_to_last_scattering = w_values;
  vector<double> kernel_values_sample_1_extended_to_last_scattering = kernel_values_sample_1;
  vector<double> kernel_values_sample_2_extended_to_last_scattering = kernel_values_sample_2;
  dw = min(constants::maximal_dw, w_last_scattering - w_values_extended_to_last_scattering[n_time]);
  w = w_values_extended_to_last_scattering[n_time] + dw;
  while(w < w_last_scattering){
    w_values_extended_to_last_scattering.push_back(w);
    kernel_values_sample_1_extended_to_last_scattering.push_back(0.0);
    kernel_values_sample_2_extended_to_last_scattering.push_back(0.0);
    n_time = w_values_extended_to_last_scattering.size()-1;
    dw = min(constants::maximal_dw, w_last_scattering - w);
    w = w + dw;
  }
  w_values_extended_to_last_scattering.push_back(w);
  kernel_values_sample_1_extended_to_last_scattering.push_back(0.0);
  kernel_values_sample_2_extended_to_last_scattering.push_back(0.0);
  n_time = w_values_extended_to_last_scattering.size()-1;
  
  vector<double> w_values_bin_center(n_time, 0.0);
  for(int i = 0; i < n_time; i++){
    w_values_bin_center[i] = 0.5*(w_values_extended_to_last_scattering[i+1]+w_values_extended_to_last_scattering[i]);
  }
  
  vector<double> CMB_lensing_kernel(n_time, 0.0);
  vector<double> CMB_lensing_kernel_overlap_sample_1(n_time, 0.0);
  vector<double> CMB_lensing_kernel_overlap_sample_2(n_time, 0.0);
  vector<double> CMB_lensing_kernel_nonoverlap(n_time, 0.0);
  double w_min = w_values_extended_to_last_scattering[0];
  double w_max = w_values_extended_to_last_scattering[n_time];
  double norm = 0.0;
  for(int i = 0; i < n_time; i++){
    w = w_values_bin_center[i];
    a = pointer_to_universe->a_at_eta(eta_0-w);
    CMB_lensing_kernel[i] = 1.5*pointer_to_universe->return_Omega_m()*w*(w_last_scattering-w)/w_last_scattering/a;
    if(kernel_values_sample_1_extended_to_last_scattering[i] > 0.0){
      CMB_lensing_kernel_overlap_sample_1[i] = CMB_lensing_kernel[i];
    }
    if(kernel_values_sample_2_extended_to_last_scattering[i] > 0.0){
      CMB_lensing_kernel_overlap_sample_2[i] = CMB_lensing_kernel[i];
    }
    if(kernel_values_sample_1_extended_to_last_scattering[i] == 0.0 && kernel_values_sample_2_extended_to_last_scattering[i] == 0.0){
      CMB_lensing_kernel_nonoverlap[i] = CMB_lensing_kernel[i];
    }
      
  }
  
  
  vector<vector<double> > joint_kernel_values(6, vector<double>(n_time, 0.0));
  
  joint_kernel_values[0] = kernel_values_sample_1_extended_to_last_scattering;
  joint_kernel_values[1] = kernel_values_sample_2_extended_to_last_scattering;
  joint_kernel_values[2] = CMB_lensing_kernel;
  joint_kernel_values[3] = CMB_lensing_kernel_overlap_sample_1;
  joint_kernel_values[4] = CMB_lensing_kernel_overlap_sample_2;
  joint_kernel_values[5] = CMB_lensing_kernel_nonoverlap;
  
  
  vector<vector<vector<double> > > C_ells = pointer_to_universe->return_LOS_integrated_C_ells(l_max, w_values_extended_to_last_scattering, joint_kernel_values);
  
  vector<vector<double> > second_moments = pointer_to_universe->return_LOS_integrated_2nd_moments(theta_in_arcmin*constants::arcmin, 0.0, 1.0, w_values_extended_to_last_scattering, joint_kernel_values);
  
  vector<vector<vector<double> > > third_moments = pointer_to_universe->return_LOS_integrated_3rd_moments(theta_in_arcmin*constants::arcmin, 0.0, 1.0, w_values_extended_to_last_scattering, joint_kernel_values);
  
  int N_ell = C_ells[0][0].size();
  
  /***** FLASK configs for sample 1: *****/
  double lambda_m_sample_1 = lognormal_tools::get_delta0(second_moments[0][0], third_moments[0][0][0]);
  double lambda_k_overlap_sample_1_correlated = lognormal_tools::get_kappa0(lambda_m_sample_1, second_moments[0][0], second_moments[0][2], third_moments[0][0][2]);
  
  double var_kappa_overlap_sample_1 = second_moments[3][3];
  double var_kappa_overlap_sample_1_correlated = lognormal_tools::get_var_kappa(lambda_m_sample_1, lambda_k_overlap_sample_1_correlated, second_moments[0][2], third_moments[0][2][2]);
  double var_kappa_overlap_sample_1_uncorrelated = var_kappa_overlap_sample_1 - var_kappa_overlap_sample_1_correlated;
  
  double skewness_kappa_overlap_sample_1 = third_moments[3][3][3];
  double skewness_kappa_overlap_sample_1_correlated = lognormal_tools::get_skew_from_delta_0(lambda_k_overlap_sample_1_correlated, var_kappa_overlap_sample_1);
  double skewness_kappa_overlap_sample_1_uncorrelated = skewness_kappa_overlap_sample_1 - skewness_kappa_overlap_sample_1_correlated;
  /***************************************/
  
  /***** FLASK configs for sample 2: *****/
  double lambda_m_sample_2 = lognormal_tools::get_delta0(second_moments[1][1], third_moments[1][1][1]);
  double lambda_k_overlap_sample_2_correlated = lognormal_tools::get_kappa0(lambda_m_sample_2, second_moments[1][1], second_moments[1][2], third_moments[1][1][2]);
  
  double var_kappa_overlap_sample_2 = second_moments[4][4];
  double var_kappa_overlap_sample_2_correlated = lognormal_tools::get_var_kappa(lambda_m_sample_2, lambda_k_overlap_sample_2_correlated, second_moments[1][2], third_moments[1][2][2]);
  double var_kappa_overlap_sample_2_uncorrelated = var_kappa_overlap_sample_2 - var_kappa_overlap_sample_2_correlated;
  
  double skewness_kappa_overlap_sample_2 = third_moments[4][4][4];
  double skewness_kappa_overlap_sample_2_correlated = lognormal_tools::get_skew_from_delta_0(lambda_k_overlap_sample_2_correlated, var_kappa_overlap_sample_2);
  double skewness_kappa_overlap_sample_2_uncorrelated = skewness_kappa_overlap_sample_2 - skewness_kappa_overlap_sample_2_correlated;
  /***************************************/
  
  double var_kappa_uncorrelated = second_moments[5][5] + var_kappa_overlap_sample_1_uncorrelated + var_kappa_overlap_sample_2_uncorrelated;
  double skewness_uncorrelated = third_moments[5][5][5] + skewness_kappa_overlap_sample_1_uncorrelated + skewness_kappa_overlap_sample_2_uncorrelated;
  double lambda_k_uncorrelated = lognormal_tools::get_delta0(var_kappa_uncorrelated, skewness_uncorrelated);
  
  
  FILE *F = fopen("FLASK_config_Cells", "w");
  fclose(F);
  fstream out;
  out.open("FLASK_config_Cells");
  out << scientific << setprecision(5);
  
  out << "# variances: ";
  out << "var_g1 = " << second_moments[0][0]*bias_sample_1*bias_sample_1;
  out << " ; var_g2 = " << second_moments[1][1]*bias_sample_2*bias_sample_2;
  out << " ;  cov_g1k = " << second_moments[0][2]*bias_sample_1*r_sample_1;
  out << " ;  cov_g2k = " << second_moments[1][2]*bias_sample_2*r_sample_2;
  out << " ;  var_k_correlated_with_g1 = " << var_kappa_overlap_sample_1_correlated;
  out << " ;  var_k_correlated_with_g2 = " << var_kappa_overlap_sample_2_correlated;
  out << " ;  var_k_uncorrelated = " << var_kappa_uncorrelated << "\n";
  
  out << "# shift params: ";
  out << "lambda_g1 = " << lambda_m_sample_1*bias_sample_1;
  out << " ; lambda_g2 = " << lambda_m_sample_2*bias_sample_2;
  out << " ; lambda_k_correlated_with_g1 = " << lambda_k_overlap_sample_1_correlated;
  out << " ; lambda_k_correlated_with_g2 = " << lambda_k_overlap_sample_2_correlated;
  out << " ; lambda_k_uncorrelated = " << lambda_k_uncorrelated << "\n";
  
  out << "# ell          C_g1g1(ell)          C_g1k(ell)          C_kk_correlated_with_g1(ell)          C_g2g2(ell)          C_g2k(ell)          C_kk_correlated_with_g2(ell)          C_kk_uncorrelated(ell)\n";
  
  for(int l = 0; l < N_ell; l++){
    out << l << setw(15);
    out << C_ells[0][0][l]*bias_sample_1*bias_sample_1 << setw(15);
    out << C_ells[0][2][l]*bias_sample_1*r_sample_1 << setw(15);
    out << C_ells[3][3][l]*var_kappa_overlap_sample_1_correlated/var_kappa_overlap_sample_1 << setw(15);
    out << C_ells[1][1][l]*bias_sample_2*bias_sample_2 << setw(15);
    out << C_ells[1][2][l]*bias_sample_2*r_sample_2 << setw(15);
    out << C_ells[4][4][l]*var_kappa_overlap_sample_2_correlated/var_kappa_overlap_sample_2 << setw(15);
    out << C_ells[5][5][l] + C_ells[3][3][l]*(1.0-var_kappa_overlap_sample_1_correlated/var_kappa_overlap_sample_1) + C_ells[4][4][l]*(1.0-var_kappa_overlap_sample_2_correlated/var_kappa_overlap_sample_2) << '\n';
  }
  
  out.close();
  
  
}



extern "C" void return_CMB_kappa_C_ells(int l_max, double bias_sample_1, double r_sample_1, int index_of_galaxy_sample){
  
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->pointer_to_universe();
  
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  
  
  int n_time = w_values.size()-1;
  double a, eta_0, w, z, dw, w_last_scattering;
  eta_0 = pointer_to_universe->eta_at_a(1.0);
  w_last_scattering = eta_0-pointer_to_universe->eta_at_a(1.0/(1.0+constants::z_last_scattering));
  
  // need to extend projection kernel to surface of last scattering:
  vector<double> w_values_extended_to_last_scattering = w_values;
  vector<double> kernel_values_extended_to_last_scattering = kernel_values;
  dw = min(constants::maximal_dw, w_last_scattering - w_values_extended_to_last_scattering[n_time]);
  w = w_values_extended_to_last_scattering[n_time] + dw;
  z = 1.0/pointer_to_universe->a_at_eta(eta_0 - w)-1.0;
  while(w < w_last_scattering && z <= 2.35){ // 2.35 is for Buzzard
    w_values_extended_to_last_scattering.push_back(w);
    kernel_values_extended_to_last_scattering.push_back(0.0);
    n_time = w_values_extended_to_last_scattering.size()-1;
    dw = min(constants::maximal_dw, w_last_scattering - w);
    w = w + dw;
    z = 1.0/pointer_to_universe->a_at_eta(eta_0 - w)-1.0;
  }
  w_values_extended_to_last_scattering.push_back(w);
  kernel_values_extended_to_last_scattering.push_back(0.0);
  n_time = w_values_extended_to_last_scattering.size()-1;
  
  vector<double> w_values_bin_center(n_time, 0.0);
  for(int i = 0; i < n_time; i++){
    w_values_bin_center[i] = 0.5*(w_values_extended_to_last_scattering[i+1]+w_values_extended_to_last_scattering[i]);
  }
  
  vector<double> CMB_lensing_kernel(n_time, 0.0);
  double w_min = w_values_extended_to_last_scattering[0];
  double w_max = w_values_extended_to_last_scattering[n_time];
  double norm = 0.0;
  for(int i = 0; i < n_time; i++){
    w = w_values_bin_center[i];
    a = pointer_to_universe->a_at_eta(eta_0-w);
    CMB_lensing_kernel[i] = 1.5*pointer_to_universe->return_Omega_m()*w*(w_last_scattering-w)/w_last_scattering/a;
  }
  
  
  vector<vector<double> > joint_kernel_values(2, vector<double>(n_time, 0.0));
  
  joint_kernel_values[0] = kernel_values_extended_to_last_scattering;
  joint_kernel_values[1] = CMB_lensing_kernel;
  
  
  vector<vector<vector<double> > > C_ells = pointer_to_universe->return_LOS_integrated_C_ells(l_max, w_values_extended_to_last_scattering, joint_kernel_values);
  int N_ell = C_ells[0][0].size();
  
  FILE *F = fopen("C_ells.dat", "w");
  fclose(F);
  fstream out;
  out.open("C_ells.dat");
  out << scientific << setprecision(5);
  out << "# ell          C_g1g1(ell)          C_g1k(ell)          C_kk(ell)\n";
  
  for(int l = 0; l < N_ell; l++){
    out << l << setw(15);
    out << C_ells[0][0][l]*bias_sample_1*bias_sample_1 << setw(15);
    out << C_ells[0][1][l]*bias_sample_1*r_sample_1 << setw(15);
    out << C_ells[1][1][l] << "\n";
  }
  
  out.close();
  
  
}



extern "C" void configure_FLASK_for_delta_g_and_kappa(int l_max, double theta_in_arcmin, double bias, double r, int index_of_lens_sample, int index_of_source_sample, const char *n_of_z_file){
  
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_lens_sample]->pointer_to_universe();

  
  vector<double> w_values;
  vector<double> w_values_1;
  vector<double> w_values_2;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  vector<double> dummy;
  
  global_universes.projected_galaxy_samples[index_of_lens_sample]->return_LOS_data(&w_values_1, &kernel_values, &dummy);
  global_universes.projected_galaxy_samples[index_of_source_sample]->return_LOS_data(&w_values_2, &dummy, &lensing_kernel_values);
  
  vector<vector<double> > rebinned_data = return_joint_binning(w_values_1, kernel_values, w_values_2, lensing_kernel_values);
  
  w_values = rebinned_data[0];
  kernel_values = rebinned_data[1];
  lensing_kernel_values = rebinned_data[2];
  int n_time = w_values.size()-1;
    
  vector<double> w_values_bin_center(n_time, 0.0);
  for(int i = 0; i < n_time; i++){
    w_values_bin_center[i] = 0.5*(w_values[i+1]+w_values[i]);
  }
  
  
  vector<double> lensing_kernel_overlap(n_time, 0.0);
  vector<double> lensing_kernel_nonoverlap(n_time, 0.0);
  for(int i = 0; i < n_time; i++){
    if(kernel_values[i] > 0.0){
      lensing_kernel_overlap[i] = lensing_kernel_values[i];
    }
    else{
      lensing_kernel_nonoverlap[i] = lensing_kernel_values[i];
    }
  }
  
  
  vector<vector<double> > joint_kernel_values(4, vector<double>(n_time, 0.0));
  
  joint_kernel_values[0] = kernel_values;
  joint_kernel_values[1] = lensing_kernel_values;
  joint_kernel_values[2] = lensing_kernel_overlap;
  joint_kernel_values[3] = lensing_kernel_nonoverlap;
  
  
  vector<vector<vector<double> > > C_ells = pointer_to_universe->return_LOS_integrated_C_ells(l_max, w_values, joint_kernel_values);
  
  vector<vector<double> > second_moments = pointer_to_universe->return_LOS_integrated_2nd_moments(theta_in_arcmin*constants::arcmin, 0.0, 1.0, w_values, joint_kernel_values);
  
  vector<vector<vector<double> > > third_moments = pointer_to_universe->return_LOS_integrated_3rd_moments(theta_in_arcmin*constants::arcmin, 0.0, 1.0, w_values, joint_kernel_values);
  
  int N_ell = C_ells[0][0].size();
  
  
  double lambda_m = lognormal_tools::get_delta0(second_moments[0][0], third_moments[0][0][0]);
  double lambda_k_overlap_correlated = lognormal_tools::get_kappa0(lambda_m, second_moments[0][0], second_moments[0][1], third_moments[0][0][1]);
  
  double var_kappa_overlap = second_moments[2][2];
  double var_kappa_overlap_correlated = lognormal_tools::get_var_kappa(lambda_m, lambda_k_overlap_correlated, second_moments[0][1], third_moments[0][1][1]);
  double var_kappa_overlap_uncorrelated = var_kappa_overlap - var_kappa_overlap_correlated;

  double skewness_kappa_overlap = third_moments[2][2][2];
  double skewness_kappa_overlap_correlated = lognormal_tools::get_skew_from_delta_0(lambda_k_overlap_correlated, var_kappa_overlap_correlated);
  double skewness_kappa_overlap_uncorrelated = skewness_kappa_overlap - skewness_kappa_overlap_correlated;
  
  double var_kappa_uncorrelated = second_moments[3][3] + var_kappa_overlap_uncorrelated;
  double skewness_uncorrelated = third_moments[3][3][3] + skewness_kappa_overlap_uncorrelated;
  double lambda_k_uncorrelated = lognormal_tools::get_delta0(var_kappa_uncorrelated, skewness_uncorrelated);
  
  
  FILE *F = fopen(n_of_z_file, "w");
  fclose(F);
  fstream out;
  out.open(string(n_of_z_file));
  out << scientific << setprecision(5);
  
  out << "# variances: var_g = ";
  out << second_moments[0][0]*bias*bias << " ;  cov_gk = ";
  out << second_moments[0][1]*bias*r << " ;  var_k_correlated = ";
  out << var_kappa_overlap_correlated << " ;  var_k_uncorrelated = ";
  out << var_kappa_uncorrelated << "\n";
  
  out << "# shift params: lambda_g = ";
  out << lambda_m*bias << " ;  lambda_k_correlated = ";
  out << lambda_k_overlap_correlated << " ;  lambda_k_uncorrelated = ";
  out << lambda_k_uncorrelated << "\n";
  
  out << "# ell          C_gg(ell)          C_gk(ell)          C_kk_correlated(ell)          C_kk_uncorrelated(ell)\n";
  
  for(int l = 0; l < N_ell; l++){
    out << l << setw(15);
    out << C_ells[0][0][l]*bias*bias << setw(15);
    out << C_ells[0][1][l]*bias*r << setw(15);
    out << C_ells[2][2][l]*var_kappa_overlap_correlated/var_kappa_overlap << setw(15);
    out << C_ells[3][3][l] + C_ells[2][2][l]*(1.0-var_kappa_overlap_correlated/var_kappa_overlap) << '\n';
  }
  
  out.close();
  
  
}



extern "C" void return_2nd_moments_2D(double* second_moments, double theta_in_arcmin, int index_of_lens_sample, int index_of_source_sample){
  
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_lens_sample]->pointer_to_universe();
  
  vector<double> w_values;
  vector<double> w_values_1;
  vector<double> w_values_2;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  vector<double> dummy;
  
  global_universes.projected_galaxy_samples[index_of_lens_sample]->return_LOS_data(&w_values_1, &kernel_values, &dummy);
  global_universes.projected_galaxy_samples[index_of_source_sample]->return_LOS_data(&w_values_2, &dummy, &lensing_kernel_values);
  
  vector<vector<double> > rebinned_data = return_joint_binning(w_values_1, kernel_values, w_values_2, lensing_kernel_values);
  
  w_values = rebinned_data[0];
  kernel_values = rebinned_data[1];
  lensing_kernel_values = rebinned_data[2];
  int n_time = w_values.size()-1;
  
  
  vector<vector<double> > joint_kernel_values(2, vector<double>(n_time, 0.0));
  
  joint_kernel_values[0] = kernel_values;
  joint_kernel_values[1] = lensing_kernel_values;
  
  vector<vector<double> > second_moments_array = pointer_to_universe->return_LOS_integrated_2nd_moments(theta_in_arcmin*constants::arcmin, 0.0, 1.0, w_values, joint_kernel_values);
  
  second_moments[0] = second_moments_array[0][0];
  second_moments[1] = second_moments_array[1][1];
  second_moments[2] = second_moments_array[0][1];
  
}





extern "C" void skewness_of_delta_L_in_cylindrical_aperture(double *skewness, double *dskewness_dR, double Omega_m, double Omega_b, double sigma_8, double n_s, double h_100, double f_NL, double R_in_Mpc_over_h, double L_in_Mpc_over_h, int PNG_modus){
  
  double a_initial = 0.000025;
  double a_final = 1.0;
  double Omega_r = 0.0; // assuming radiation content is negligible
  double Omega_L = 1.0 - Omega_m; // assuming spatially flat Universe
  double w0 = -1.0; // assuming dark energy is cosmological constant
  double w1 = 0.0; // assuming dark energy is cosmological constant
  
  cosmological_model cosmo = set_cosmological_model(Omega_m, Omega_b, Omega_r, Omega_L, sigma_8, n_s, h_100, w0, w1);
  
  FlatInhomogeneousUniverseLCDM our_universe(cosmo, a_initial, a_final, 1);
  
  our_universe.compute_cylinder_skewnesses_for_unit_L_and_unit_fNL(PNG_modus, R_in_Mpc_over_h, skewness, dskewness_dR);
  
  (*skewness) *= f_NL;
  (*dskewness_dR) *= f_NL;
  
  // the function compute_cylinder_skewnesses_for_unit_L_and_unit_fNL only returns
  // skewness for cylinder of length c/H_0 ==> re-scale to desired length
  // assuming L >> R .
  (*skewness) *= pow(constants::c_over_e5/L_in_Mpc_over_h, 2.0);
  (*dskewness_dR) *= pow(constants::c_over_e5/L_in_Mpc_over_h, 2.0);
  
}



extern "C" double return_R_in_Mpc_over_h_from_angular_scale_and_z(int index_of_universe, double z, double theta_in_arcmin){
  double e_0 = global_universes.universes[index_of_universe]->eta_at_a(1.0);
  double w = e_0 - global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z));
  return w*theta_in_arcmin*constants::arcmin*constants::c_over_e5;
}

extern "C" double return_angular_scale_in_arcmin_from_physical_scale_and_z(int index_of_universe, double z, double R_in_Mpc_over_h){
  double e_0 = global_universes.universes[index_of_universe]->eta_at_a(1.0);
  double w = e_0 - global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z));
  return (R_in_Mpc_over_h/constants::c_over_e5)/w/constants::arcmin;
}


extern "C" double return_Mass_within_R_in_Mpc_over_h(int index_of_universe, double R_in_Mpc_over_h){
  return global_universes.universes[index_of_universe]->Mass_within_R_in_Mpc_over_h(R_in_Mpc_over_h);
}


extern "C" double return_R_in_Mpc_over_h_enclosing_M(int index_of_universe, double M_in_Msol_over_h){
  return global_universes.universes[index_of_universe]->R_in_Mpc_over_h_enclosing_M(M_in_Msol_over_h);
}


extern "C" double return_delta_crit(int index_of_universe, double z, double Delta){
  double a = 1.0/(1.0+z);
  return global_universes.universes[index_of_universe]->delta_crit(a, Delta);
}




/*
 * PaNaMO hack week
 * 
 */




extern "C" void test_interface(double* delta_values, double* kappa_values, double* PDF, double theta_in_arcmin, double var_kappa_noise, int index_of_lens_sample, int index_of_source_sample){
  
  cout << "Universe creation started\n";
  FlatInhomogeneousUniverseLCDM* pointer_to_universe = global_universes.projected_galaxy_samples[index_of_lens_sample]->pointer_to_universe();
  
  cout << "Universe created\n";
  
  vector<double> w_values;
  vector<double> w_values_1;
  vector<double> w_values_2;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  vector<double> dummy;
  
  global_universes.projected_galaxy_samples[index_of_lens_sample]->return_LOS_data(&w_values_1, &kernel_values, &dummy);
  global_universes.projected_galaxy_samples[index_of_source_sample]->return_LOS_data(&w_values_2, &dummy, &lensing_kernel_values);
  
  vector<vector<double> > rebinned_data = return_joint_binning(w_values_1, kernel_values, w_values_2, lensing_kernel_values);
  
  w_values = rebinned_data[0];
  kernel_values = rebinned_data[1];
  lensing_kernel_values = rebinned_data[2];
  
  double f_NL = 0.0;
  double var_NL_rescale = 1.0;
  cout << "start of return_LOS_joint_PDF_deltaG_kappa\n";
  vector<vector<double> > PDF_data = pointer_to_universe->return_LOS_joint_PDF_deltaG_kappa(theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale, var_kappa_noise, w_values, kernel_values, lensing_kernel_values);
  
  int n_delta = constants::N_delta_values_for_PDFs;
  int n_kappa = constants::N_delta_values_for_PDFs;
  
  for(int id = 0; id < n_delta; id++){
    delta_values[id] = PDF_data[id+1][0];
  }
  
  for(int ik = 0; ik < n_kappa; ik++){
    kappa_values[ik] = PDF_data[0][ik+1];
  }
  
  for(int id = 0; id < n_delta; id++){
    for(int ik = 0; ik < n_kappa; ik++){
      PDF[id*n_kappa + ik] = PDF_data[id+1][ik+1];
    }
  }
  
  
}


extern "C" void set_cosmic_power_grid(double* k_vals, double* z_vals, double* P_L_kz_vals, double* P_NL_kz_vals, int n_k, int n_z, int index_of_universe){
  
  vector<double> k_values(n_k, 0.0);
  vector<double> z_values(n_z, 0.0);
  vector<vector<double> > P_L_kz_values(n_k, vector<double>(n_z, 0.0));
  vector<vector<double> > P_NL_kz_values(n_k, vector<double>(n_z, 0.0));
  
  for(int i = 0; i < n_k; i++){
    k_values[i] = k_vals[i]*constants::c_over_e5;
  }
  for(int i = 0; i < n_z; i++){
    z_values[i] = z_vals[i];
  }
  for(int i = 0; i < n_k; i++){
    for(int j = 0; j < n_z; j++){
      P_L_kz_values[i][j] = P_L_kz_vals[i*n_z + j]/pow(constants::c_over_e5,3);
      P_NL_kz_values[i][j] = P_NL_kz_vals[i*n_z + j]/pow(constants::c_over_e5,3);
    }
  }
  
  global_universes.universes[index_of_universe]->set_power_grid(&k_values, &z_values, &P_L_kz_values, &P_NL_kz_values);
  
}



