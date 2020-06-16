
#include "constants.h"
#include "lognormal_tools.h"
#include "io_utils.h"
#include "error_handling.h"
#include "interpolators.h"
#include "FlatInhomogeneousUniverseLCDM.h"
#include "GalaxySample3D.h"
#include "ProjectedGalaxySample.h"


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


extern "C" void add_3D_galaxy_sample(int index_of_universe, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1){
  global_universes.galaxy_samples_3D.push_back(new GalaxySample3D(global_universes.universes[index_of_universe], z, density_in_Mpc_over_h_cubed, b1, b2, a0, a1));
}

extern "C" int return_N_max_3D(int index_of_galaxy_sample, double R_in_Mpc_over_h, double var_NL_rescale){
  return global_universes.galaxy_samples_3D[index_of_galaxy_sample]->return_N_max_in_3D_tophat(R_in_Mpc_over_h, var_NL_rescale);
}

extern "C" void change_parameters_of_3D_galaxy_sample(int index_of_galaxy_sample, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1){
  global_universes.galaxy_samples_3D[index_of_galaxy_sample]->set_parameters_3D(z, density_in_Mpc_over_h_cubed, b1, b2, a0, a1);
}

extern "C" double change_b2_to_minimise_negative_densities_3D(int index_of_galaxy_sample, double R_in_Mpc_over_h, double var_NL_rescale){
  return global_universes.galaxy_samples_3D[index_of_galaxy_sample]->set_b2_to_minimise_negative_densities_in_3D_tophat(R_in_Mpc_over_h, var_NL_rescale);
}

extern "C" void update_3D_bias_model_from_br_parametrisation(int index_of_galaxy_sample, double b_tilde, double r, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale){
  global_universes.galaxy_samples_3D[index_of_galaxy_sample]->set_3D_bias_model_from_br_parametrisation(b_tilde, r, R_in_Mpc_over_h, f_NL, var_NL_rescale);
}

extern "C" void return_CiC_PDF_3D(double* P_of_N, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  vector<double> P_of_N_vector = global_universes.galaxy_samples_3D[index_of_galaxy_sample]->return_CiC_PDF_in_3D_tophat(R_in_Mpc_over_h, f_NL, var_NL_rescale);
  
  int number_of_Ns = P_of_N_vector.size();
  
  for(int i = 0; i < number_of_Ns; i++){
    P_of_N[i] = P_of_N_vector[i];
  }
  
}


extern "C" void add_projected_galaxy_sample(int index_of_universe, const char *n_of_z_file, double density_in_arcmin_squared, double b1, double b2, double a0, double a1){
  global_universes.projected_galaxy_samples.push_back(new ProjectedGalaxySample(global_universes.universes[index_of_universe], density_in_arcmin_squared, b1, b2, a0, a1, string(n_of_z_file)));
}

extern "C" int return_N_max_projected(int index_of_galaxy_sample, double theta, double var_NL_rescale){
  return global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_N_max_in_angular_tophat(theta, var_NL_rescale);
}

extern "C" void change_parameters_of_projected_galaxy_sample(int index_of_galaxy_sample, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, const char *n_of_z_file){
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->set_parameters_projected(density_in_arcmin_squared, b1, b2, a0, a1, string(n_of_z_file));
}

extern "C" double change_b2_to_minimise_negative_densities_projected(int index_of_galaxy_sample, double theta, double var_NL_rescale){
  return global_universes.projected_galaxy_samples[index_of_galaxy_sample]->set_b2_to_minimise_negative_densities_in_angular_tophat(theta, var_NL_rescale);
}

extern "C" void return_CiC_PDF_projected(double* P_of_N, double theta, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  vector<double> P_of_N_vector = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_CiC_PDF_in_angular_tophat(theta, f_NL, var_NL_rescale);
  
  int number_of_Ns = P_of_N_vector.size();
  
  for(int i = 0; i < number_of_Ns; i++){
    P_of_N[i] = P_of_N_vector[i];
  }
  
}

extern "C" void return_CiC_saddle_point_PDF_projected(double* P_of_N, double theta, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  vector<double> P_of_N_vector = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_CiC_saddle_point_PDF_in_angular_tophat(theta, f_NL, var_NL_rescale);
  
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

extern "C" void return_power_spectra(double* k_values, double* P_L, double* P_halofit, int index_of_universe, double z, double R){
  
  double eta = global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z));
  
  vector<vector<double> > power_spectra = global_universes.universes[index_of_universe]->return_power_spectra(eta, R);
  
  for(int i = 0; i < constants::number_of_k; i++){
    k_values[i] = power_spectra[0][i];
    P_L[i]   = power_spectra[1][i];
    P_halofit[i]   = power_spectra[2][i];
  }
  
}

extern "C" void return_CGF(double* delta_L, double* delta_NL, double* lambda, double* phi, double* lambda_Gauss, double* phi_Gauss, double* variances, double* skewnesses, double* R_L, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > phi_data = global_universes.universes[index_of_universe]->compute_phi_of_lambda_3D(z, R_in_comoving_Mpc/constants::c_over_e5, f_NL, var_NL_rescale);

  for(int i = 0; i < phi_data[0].size(); i++){
    delta_L[i] = phi_data[0][i];
    delta_NL[i] = phi_data[1][i];
    lambda[i] = phi_data[2][i];
    phi[i] = phi_data[3][i];
    lambda_Gauss[i] = phi_data[4][i];
    phi_Gauss[i] = phi_data[5][i];
    variances[i] = phi_data[6][i];
    skewnesses[i] = phi_data[7][i];
    R_L[i] = phi_data[8][i];
  }
  
}

extern "C" void return_PDF(double* delta_values, double* PDF, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > PDF_data = global_universes.universes[index_of_universe]->compute_PDF_3D(z, R_in_comoving_Mpc, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    delta_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
  }
  
}

extern "C" int return_Ndelta(){
  return constants::N_delta_values_for_PDFs;
}

extern "C" int return_Nlambda(int index_of_universe){
  return global_universes.universes[index_of_universe]->return_N_of_lambda();
}

extern "C" double return_var_NL(double z, double R_in_Mpc_over_h, int index_of_universe){
  // ISSUE: var_NL_rescale non implemented here
  return global_universes.universes[index_of_universe]->return_non_linear_variance(z, R_in_Mpc_over_h, 1.0);
}


extern "C" double return_var_L(double z, double R_in_Mpc_over_h, int index_of_universe){
  return global_universes.universes[index_of_universe]->return_linear_variance(z, R_in_Mpc_over_h);
}

extern "C" void compute_projected_PDF_incl_CMB_kappa(int *Nd, int *Nk, double theta_in_arcmin, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  vector<vector<double> > *d_grid = &global_universes.delta_grid;
  vector<vector<double> > *k_grid = &global_universes.kappa_grid;
  vector<vector<double> > *p_grid = &global_universes.PDF_grid;
  vector<double> w_values;
  vector<double> kernel_values;
  vector<double> lensing_kernel_values;
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_LOS_data(&w_values, &kernel_values, &lensing_kernel_values);
  global_universes.projected_galaxy_samples[index_of_galaxy_sample]->pointer_to_universe()->compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point(w_values, kernel_values, theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale, d_grid, k_grid, p_grid);
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

extern "C" void return_convergence_PDF_from_source_sample(double* kappa_values, double* PDF, double theta_in_arcmin, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  vector<vector<double> > PDF_data = global_universes.projected_galaxy_samples[index_of_galaxy_sample]->return_kappa_PDF_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale);
  
  for(int d = 0; d < PDF_data[0].size(); d++){
    kappa_values[d] = PDF_data[0][d];
    PDF[d] = PDF_data[1][d];
  }
  
}


