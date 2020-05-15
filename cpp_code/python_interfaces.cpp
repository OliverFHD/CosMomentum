
#include "constants.h"
#include "io_utils.h"
#include "error_handling.h"
#include "interpolators.h"
#include "Universe.cpp"
#include "Matter.cpp"
#include "GalaxySample.cpp"
#include "GalaxySample3D.cpp"
#include "ProjectedGalaxySample.cpp"


class Global_Universes {

  public:
    Global_Universes(){}
    
    cosmological_model cosmo;
    
    vector<Universe*> universes;
    vector<Matter*> matter_contents;
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
  
  int does_universe_expand = 1; // 1 = yes; 0 = no
  if(a_final < a_initial) does_universe_expand = 0;
  
  global_universes.universes.push_back(new Universe(cosmo, a_initial, a_final, does_universe_expand));
  
}


extern "C" void initialise_new_Universe_with_Matter_content(double a_initial, double a_final, double Omega_m, double Omega_b, double Omega_r, double Omega_L, double sigma_8, double n_s, double h_100, double w0, double w1){
  
  cosmological_model cosmo = set_cosmological_model(Omega_m, Omega_b, Omega_r, Omega_L, sigma_8, n_s, h_100, w0, w1);
  global_universes.cosmo = cosmo;
  
  int does_universe_expand = 1; // 1 = yes; 0 = no
  if(a_final < a_initial) does_universe_expand = 0;
    
  global_universes.universes.push_back(new Universe(cosmo, a_initial, a_final, does_universe_expand));
  int N_Universe = global_universes.universes.size();
  global_universes.matter_contents.push_back(new Matter(global_universes.universes[N_Universe-1]));
  //global_universes.matter_contents.push_back(new Matter(global_universes.universes[N_Universe-1], "../Density_Split_Statistics/code/Coras_PDFs/input_transfer_manera.dat"));
  int N_Matter = global_universes.matter_contents.size();
  
  error_handling::test_value_of_int(N_Universe, N_Matter, error_handling::EQUAL, "Number of Universes = " + to_string(N_Universe) + " does not match Number of Matter contents = " + to_string(N_Matter) + " .");
  
}


extern "C" void add_3D_galaxy_sample(int index_of_universe, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1){
  global_universes.galaxy_samples_3D.push_back(new GalaxySample3D(global_universes.matter_contents[index_of_universe], z, density_in_Mpc_over_h_cubed, b1, b2, a0, a1));
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

extern "C" void return_CiC_PDF_3D(double* P_of_N, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale, int index_of_galaxy_sample){
  
  vector<double> P_of_N_vector = global_universes.galaxy_samples_3D[index_of_galaxy_sample]->return_CiC_PDF_in_3D_tophat(R_in_Mpc_over_h, f_NL, var_NL_rescale);
  
  int number_of_Ns = P_of_N_vector.size();
  
  for(int i = 0; i < number_of_Ns; i++){
    P_of_N[i] = P_of_N_vector[i];
  }
  
}


extern "C" void add_projected_galaxy_sample(int index_of_universe, const char *n_of_z_file, double density_in_arcmin_squared, double b1, double b2, double a0, double a1){
  global_universes.projected_galaxy_samples.push_back(new ProjectedGalaxySample(global_universes.matter_contents[index_of_universe], density_in_arcmin_squared, b1, b2, a0, a1, string(n_of_z_file)));
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



extern "C" void set_primordial_skewness(int index_of_universe, int PNG_modus){
  cout << "Initialising primordial skewnesses from exact integrals.\n";
  cout << "NOTE: this can be very slow!";
  global_universes.matter_contents[index_of_universe]->set_sphere_skewnesses(PNG_modus);
}

extern "C" void set_primordial_skewness_from_eps3_powerlaw_approximation(int index_of_universe, int PNG_modus, double R_0_in_Mpc_over_h){
  cout << "Initialising primordial skewnesses from power law approximation\n";
  cout << "eps3(R) = skewness(R)/sigma(R)^3 = A_eps3*(R/R_0)^n_eps3 .\n";
  global_universes.matter_contents[index_of_universe]->set_sphere_skewnesses_from_eps3_powerlaw_approximation(PNG_modus, R_0_in_Mpc_over_h);
}

// skewnesses stored in "skewness_file" should be for f_NL = 1;
extern "C" void set_primordial_skewness_from_file(int index_of_universe, const char *skewness_file){
  cout << "Initialising primordial skewnesses from file\n";
  cout << skewness_file << " .\n";
  cout << "NOTE: skewnesses in this file should be for f_NL = 1";
  global_universes.matter_contents[index_of_universe]->set_sphere_skewnesses_from_file(skewness_file);
}




extern "C" void set_primordial_skewness_2D(int index_of_universe, int PNG_modus){
  cout << "Initialising primordial skewnesses from exact integrals.\n";
  cout << "NOTE: this can be very slow!";
  global_universes.matter_contents[index_of_universe]->set_cylinder_skewnesses(PNG_modus);
}

extern "C" void set_primordial_skewness_from_eps3_powerlaw_approximation_2D(int index_of_universe, int PNG_modus, double R_0_in_Mpc_over_h){
  cout << "Initialising primordial skewnesses from power law approximation\n";
  cout << "eps3(R) = skewness(R)/sigma(R)^3 = A_eps3*(R/R_0)^n_eps3 .\n";
  global_universes.matter_contents[index_of_universe]->set_cylinder_skewnesses_from_eps3_powerlaw_approximation(PNG_modus, R_0_in_Mpc_over_h);
}

extern "C" void clear_universes(){
  
  for(int i = 0; i < global_universes.universes.size();i++) delete global_universes.universes[i];
  for(int i = 0; i < global_universes.matter_contents.size();i++) delete global_universes.matter_contents[i];
  for(int i = 0; i < global_universes.galaxy_samples_3D.size();i++) delete global_universes.galaxy_samples_3D[i];
  for(int i = 0; i < global_universes.projected_galaxy_samples.size();i++) delete global_universes.projected_galaxy_samples[i];
  
  global_universes.universes.clear();
  global_universes.matter_contents.clear();
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
  
  vector<vector<double> > linear_growth_history = global_universes.matter_contents[index_of_universe]->return_linear_growth_history(conformal_time_steps);
  
  int n_column = linear_growth_history.size();
  
  for(int i = 0; i < conformal_time_steps; i++){
    eta[i] = linear_growth_history[0][i];
    a[i]   = linear_growth_history[1][i];
    D[i]   = linear_growth_history[2][i];
    D_prime[i] = linear_growth_history[3][i];
  }
  
}


extern "C" void return_flat_matter_dominated_background_expansion_to_python(double* t_phys, double* eta, double* a, double* H_conf, double* H_conf_prime, int conformal_time_steps){
  
  for(int i = 0; i < conformal_time_steps; i++) Universe::expansion_in_flat_matter_dominated_universe(a[i], 1.0, &t_phys[i], &eta[i], &H_conf[i], &H_conf_prime[i]);  
  
}

extern "C" void return_flat_radiation_dominated_background_expansion_to_python(double* t_phys, double* eta, double* a, double* H_conf, double* H_conf_prime, int conformal_time_steps){
  
  for(int i = 0; i < conformal_time_steps; i++) Universe::expansion_in_flat_radiation_dominated_universe(a[i], &t_phys[i], &eta[i], &H_conf[i], &H_conf_prime[i]);  
  
}

extern "C" void return_flat_Lambda_dominated_background_expansion_to_python(double* t_phys, double* eta, double* a, double* H_conf, double* H_conf_prime, int conformal_time_steps){
  
  for(int i = 0; i < conformal_time_steps; i++) Universe::expansion_in_flat_Lambda_dominated_universe(a[i], &t_phys[i], &eta[i], &H_conf[i], &H_conf_prime[i]);  
  
}

extern "C" int return_n_k(){
  return constants::number_of_k;
}

extern "C" void return_power_spectra(double* k_values, double* P_L, double* P_halofit, int index_of_universe, double z, double R){
  
  double eta = global_universes.universes[index_of_universe]->eta_at_a(1.0/(1.0+z));
  
  vector<vector<double> > power_spectra = global_universes.matter_contents[index_of_universe]->return_power_spectra(eta, R);
  
  for(int i = 0; i < constants::number_of_k; i++){
    k_values[i] = power_spectra[0][i];
    P_L[i]   = power_spectra[1][i];
    P_halofit[i]   = power_spectra[2][i];
  }
  
}

extern "C" void return_CGF(double* delta_L, double* delta_NL, double* lambda, double* phi, double* lambda_Gauss, double* phi_Gauss, double* variances, double* skewnesses, double* R_L, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > phi_data = global_universes.matter_contents[index_of_universe]->compute_phi_of_lambda_3D(z, R_in_comoving_Mpc/constants::c_over_e5, f_NL, var_NL_rescale);

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

extern "C" void return_PDF(double* delta_values, double* PDF, int n_delta, double delta_min, double delta_max, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  cout << "Computing phi_data:\n";
  cout.flush();
  
  vector<vector<double> > phi_data;
  
  phi_data = global_universes.matter_contents[index_of_universe]->compute_phi_of_lambda_3D(z, R_in_comoving_Mpc/constants::c_over_e5, f_NL, var_NL_rescale);
  
  /*
   * Determine critical point, where phi(lambda) splits into two branches on the complex plane. 
   * 
   */
  int n_lambda = 0;
  double lambda_c = phi_data[2][0];
  double delta_c = phi_data[1][0];
  double tau_c = phi_data[0][0]/sqrt(phi_data[6][0]);
  for(int i = 1; i < phi_data[4].size(); i++){
    if(phi_data[4][i-1] < phi_data[4][i]){
      n_lambda = i+1;
      lambda_c = phi_data[4][i];
      delta_c = phi_data[1][i];
      tau_c = phi_data[0][i]/sqrt(phi_data[6][i]);
    }
    else{
      i = 2*phi_data[4].size();
    }
  }
  
  /*
   * Extract phi_data up to the critical point.
   * 
   */
  vector<double> delta_L_values(n_lambda, 0.0);
  vector<double> delta_NL_values(n_lambda, 0.0);
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> phi_values(n_lambda, 0.0);
  vector<double> lambda_values_Gauss(n_lambda, 0.0);
  vector<double> phi_values_Gauss(n_lambda, 0.0);
  vector<double> variances(n_lambda, 0.0);
  vector<double> skewnesses(n_lambda, 0.0);
  vector<double> R_L_values(n_lambda, 0.0);
  vector<double> tau_values(n_lambda, 0.0);
  
  for(int i = 0; i < tau_values.size(); i++){
    delta_L_values[i] = phi_data[0][i];
    delta_NL_values[i] = phi_data[1][i];
    lambda_values[i] = phi_data[2][i];
    phi_values[i] = phi_data[3][i];
    lambda_values_Gauss[i] = phi_data[4][i];
    phi_values_Gauss[i] = phi_data[5][i];
    variances[i] = phi_data[6][i];
    skewnesses[i] = phi_data[7][i];
    R_L_values[i] = phi_data[8][i];
    tau_values[i] = delta_L_values[i]/sqrt(variances[i]);
  }
  
  
  /*
   * Extract phi_data with equal number and range of points left and right of tau = 0 (better for polynomial fit).
   * 
   */
  
  double tau_max = 0.9*tau_c;
  double tau_min = 0.9*tau_values[0];
  
  double R_max = interpolate_neville_aitken(tau_max, &tau_values, &R_L_values, constants::order_of_interpolation);
  
  double R_min = interpolate_neville_aitken(tau_min, &tau_values, &R_L_values, constants::order_of_interpolation);
    
  if(R_max > exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5){
    tau_max = min(tau_max, interpolate_neville_aitken(exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  }
  
  if(R_min < exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5){
    tau_min = max(tau_min, interpolate_neville_aitken(exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  }
  
  tau_c = tau_max;
    
  delta_c = interpolate_neville_aitken(tau_max, &tau_values, &delta_NL_values, constants::order_of_interpolation);
  lambda_c = interpolate_neville_aitken(tau_max, &tau_values, &lambda_values, constants::order_of_interpolation);
  
  int n_tau = 4*constants::generating_function_coeff_order + 1; // has to be odd number in order to include tau=0 exactly.
  
  vector<double> tau_for_fit(n_tau,0.0);
  vector<double> lambda_for_fit(n_tau,0.0);
  vector<double> phi_for_fit(n_tau,0.0);
  vector<double> phi_prime_for_fit(n_tau,0.0);
  
  
  double dt = -tau_min/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i] = tau_min+double(i)*dt;
  }
  tau_for_fit[n_tau/2] = 0.0;
  dt = tau_max/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i+n_tau/2+1] = double(i+1)*dt;
  }
  
  for(int i = 0; i < n_tau; i++){
    lambda_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &lambda_values, constants::order_of_interpolation);
    phi_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &phi_values, constants::order_of_interpolation);
    phi_prime_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &delta_NL_values, constants::order_of_interpolation);
  }
  
  cout << "Done.\n";
  
  /*
   * Express functions as polynomials in tau.
   * 
   */
  cout << "Computing tau coefficients:\n";
  
  int n_coeff = constants::generating_function_coeff_order;
  
  vector<double> coefficients_lambda_of_tau = return_coefficients(&tau_for_fit, &lambda_for_fit, n_coeff);
  vector<double> coefficients_lambda_of_tau_prime(coefficients_lambda_of_tau.size(), 0.0);
  for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++) coefficients_lambda_of_tau_prime[i] = coefficients_lambda_of_tau[i+1]*double(i+1);
  
  vector<double> coefficients_phi_of_tau = return_coefficients(&tau_for_fit, &phi_for_fit, n_coeff);
  vector<double> coefficients_phi_of_tau_prime = return_coefficients(&tau_for_fit, &phi_prime_for_fit, n_coeff);
  
  cout << "Done.\n";
  
  /*
   * Perform the inverse Laplace transform of phi(lambda) to compute p(delta).
   * 
   */
  double ddelta = (delta_max-delta_min)/double(n_delta-1);
  double delta;
  double tau_0, lambda_0;
  double dr;
  
  complex<double> lambda;
  complex<double> lambda_next;
  complex<double> tau, tau_next;
  complex<double> phi_prime, phi_prime_next;
  complex<double> exponent, exponent_next;
  complex<double> dlambda;
  complex<double> step, first_step;
  
  cout << "Computing PDF:\n";
  for(int i = 0; i < n_delta; i++){
    delta = delta_min + double(i)*ddelta;
    
    delta_values[i] = delta;
    PDF[i] = 0.0;
      
    if(delta < delta_c){
      tau_0 = interpolate_Newton(delta, &delta_NL_values, &tau_values, constants::order_of_interpolation);
      lambda_0 = interpolate_Newton(delta, &delta_NL_values, &lambda_values, constants::order_of_interpolation);
    }
    else{
      tau_0 = tau_c;
      lambda_0 = lambda_c;
    }
    lambda = complex<double>(lambda_0, 0.0);
    tau = complex<double>(tau_0, 0.0);
    exponent = exp(-lambda*delta + return_polnomial_value(tau, &coefficients_phi_of_tau));
    
    // sigma_r^2 \approx 1/phi''(lambda_0)
    double sigma_frac = 0.001;
    dr = sigma_frac/sqrt(interpolate_neville_aitken_derivative(lambda_0, &lambda_values, &delta_NL_values, constants::order_of_interpolation));
    dlambda = complex<double>(0.0, dr);
    int j = 0;
    
    //if(i%1 == 0){
    //  cout << delta << setw(20);
    //  cout.flush();
    //}
    do{
      lambda_next = lambda + 0.5*dlambda;
      tau_next = Newtons_method_complex(tau, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda_next = lambda + dlambda;
      tau_next = Newtons_method_complex(tau_next, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      exponent_next = exp(-lambda_next*delta + return_polnomial_value(tau_next, &coefficients_phi_of_tau));
      
      step = 0.5*dlambda*(exponent_next+exponent);
      PDF[i] += step.imag();
      
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda = lambda_next;
      tau = tau_next;
      exponent = exponent_next;
      if(j == 0){
        first_step = step;
      }
      j++;
    }while(abs(step/first_step) > 1.0e-5 || j < int(6.0/sigma_frac));
    
    PDF[i] /= constants::pi;
    
    //if(i%1 == 0){
    //  cout << PDF[i] << setw(20);
    //  cout  << '\n';
    //}
    
  }
  cout << "Done.\n";
  
}






extern "C" void return_PDF_tauEff_instead_of_tau(double* delta_values, double* PDF, int n_delta, double delta_min, double delta_max, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  cout << "Computing phi_data:\n";
  cout.flush();
  
  vector<vector<double> > phi_data;
  
  phi_data = global_universes.matter_contents[index_of_universe]->compute_phi_of_lambda_3D(z, R_in_comoving_Mpc/constants::c_over_e5, f_NL, var_NL_rescale);
  
  /*
   * Determine critical point, where phi(lambda) splits into two branches on the complex plane. 
   * 
   */
  int n_lambda = 0;
  double lambda_c = phi_data[2][0];
  double delta_c = phi_data[1][0];
  double tau_c = 0.0;
  for(int i = 1; i < phi_data[4].size(); i++){
    if(phi_data[4][i-1] < phi_data[4][i]){
      n_lambda = i+1;
      lambda_c = phi_data[4][i];
      delta_c = phi_data[1][i];
      if(phi_data[1][i]*phi_data[4][i] > phi_data[3][i])
        tau_c = sqrt(2.0*(phi_data[1][i]*phi_data[4][i] - phi_data[3][i]));
    }
    else{
      i = 2*phi_data[4].size();
    }
  }
  
  /*
   * Extract phi_data up to the critical point.
   * 
   */
  vector<double> delta_L_values(n_lambda, 0.0);
  vector<double> delta_NL_values(n_lambda, 0.0);
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> phi_values(n_lambda, 0.0);
  vector<double> lambda_values_Gauss(n_lambda, 0.0);
  vector<double> phi_values_Gauss(n_lambda, 0.0);
  vector<double> variances(n_lambda, 0.0);
  vector<double> skewnesses(n_lambda, 0.0);
  vector<double> R_L_values(n_lambda, 0.0);
  vector<double> tau_values(n_lambda, 0.0);
  
  for(int i = 0; i < tau_values.size(); i++){
    delta_L_values[i] = phi_data[0][i];
    delta_NL_values[i] = phi_data[1][i];
    lambda_values[i] = phi_data[2][i];
    phi_values[i] = phi_data[3][i];
    lambda_values_Gauss[i] = phi_data[4][i];
    phi_values_Gauss[i] = phi_data[5][i];
    variances[i] = phi_data[6][i];
    skewnesses[i] = phi_data[7][i];
    R_L_values[i] = phi_data[8][i];
    if(phi_data[1][i]*phi_data[2][i] > phi_data[3][i]){
      tau_values[i] = sqrt(2.0*(phi_data[1][i]*phi_data[2][i] - phi_data[3][i]));
      if(lambda_values[i] < 0.0) tau_values[i] *= -1.0;
    }
  }
  
  
  /*
   * Extract phi_data with equal number and range of points left and right of tau = 0 (better for polynomial fit).
   * 
   */
  
  double tau_max = 0.9*tau_c;
  double tau_min = 0.9*tau_values[0];
  
  double R_max = interpolate_neville_aitken(tau_max, &tau_values, &R_L_values, constants::order_of_interpolation);
  
  double R_min = interpolate_neville_aitken(tau_min, &tau_values, &R_L_values, constants::order_of_interpolation);
  
  cout << tau_min << "  " << tau_max << "  " << R_min << "  " << R_max << "\n";
  
  if(R_max > exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5){
    tau_max = min(tau_max, interpolate_neville_aitken(exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  }
  
  if(R_min < exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5){
    tau_min = max(tau_min, interpolate_neville_aitken(exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  }
  
  tau_c = tau_max;
  
  cout << tau_min << "  " << tau_max << "  " << R_min << "  " << R_max << "\n";
  
  delta_c = interpolate_neville_aitken(tau_max, &tau_values, &delta_NL_values, constants::order_of_interpolation);
  lambda_c = interpolate_neville_aitken(tau_max, &tau_values, &lambda_values, constants::order_of_interpolation);
  
  int n_tau = 4*constants::generating_function_coeff_order + 1; // has to be odd number in order to include tau=0 exactly.
  
  vector<double> tau_for_fit(n_tau,0.0);
  vector<double> lambda_for_fit(n_tau,0.0);
  vector<double> phi_for_fit(n_tau,0.0);
  vector<double> phi_prime_for_fit(n_tau,0.0);
  
  
  double dt = -tau_min/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i] = tau_min+double(i)*dt;
  }
  tau_for_fit[n_tau/2] = 0.0;
  dt = tau_max/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i+n_tau/2+1] = double(i+1)*dt;
  }
  
  for(int i = 0; i < n_tau; i++){
    lambda_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &lambda_values, constants::order_of_interpolation);
    phi_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &phi_values, constants::order_of_interpolation);
    phi_prime_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &delta_NL_values, constants::order_of_interpolation);
  }
  
  cout << "Done.\n";
  
  /*
   * Express functions as polynomials in tau.
   * 
   */
  cout << "Computing tau coefficients:\n";
  
  int n_coeff = constants::generating_function_coeff_order;
  
  vector<double> coefficients_lambda_of_tau = return_coefficients(&tau_for_fit, &lambda_for_fit, n_coeff);
  vector<double> coefficients_lambda_of_tau_prime(coefficients_lambda_of_tau.size(), 0.0);
  for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++) coefficients_lambda_of_tau_prime[i] = coefficients_lambda_of_tau[i+1]*double(i+1);
  
  vector<double> coefficients_phi_of_tau = return_coefficients(&tau_for_fit, &phi_for_fit, n_coeff);
  vector<double> coefficients_phi_of_tau_prime = return_coefficients(&tau_for_fit, &phi_prime_for_fit, n_coeff);
  
  cout << "Done.\n";
  
  /*
   * Perform the inverse Laplace transform of phi(lambda) to compute p(delta).
   * 
   */
  double ddelta = (delta_max-delta_min)/double(n_delta-1);
  double delta;
  double tau_0, lambda_0;
  double dr;
  
  complex<double> lambda;
  complex<double> lambda_next;
  complex<double> tau, tau_next;
  complex<double> phi_prime, phi_prime_next;
  complex<double> exponent, exponent_next;
  complex<double> dlambda;
  complex<double> step, first_step;
  
  cout << "Computing PDF:\n";
  for(int i = 0; i < n_delta; i++){
    delta = delta_min + double(i)*ddelta;
    
    delta_values[i] = delta;
    PDF[i] = 0.0;
      
    if(delta < delta_c){
      tau_0 = interpolate_Newton(delta, &delta_NL_values, &tau_values, constants::order_of_interpolation);
      lambda_0 = interpolate_Newton(delta, &delta_NL_values, &lambda_values, constants::order_of_interpolation);
    }
    else{
      tau_0 = tau_c;
      lambda_0 = lambda_c;
    }
    lambda = complex<double>(lambda_0, 0.0);
    tau = complex<double>(tau_0, 0.0);
    exponent = exp(-lambda*delta + return_polnomial_value(tau, &coefficients_phi_of_tau));
    
    // sigma_r^2 \approx 1/phi''(lambda_0)
    double sigma_frac = 0.001;
    dr = sigma_frac/sqrt(interpolate_neville_aitken_derivative(lambda_0, &lambda_values, &delta_NL_values, constants::order_of_interpolation));
    dlambda = complex<double>(0.0, dr);
    int j = 0;
    
    if(i%1 == 0){
      cout << delta << setw(20);
      cout.flush();
    }
    do{
      lambda_next = lambda + 0.5*dlambda;
      tau_next = Newtons_method_complex(tau, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda_next = lambda + dlambda;
      tau_next = Newtons_method_complex(tau_next, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      exponent_next = exp(-lambda_next*delta + return_polnomial_value(tau_next, &coefficients_phi_of_tau));
      
      step = 0.5*dlambda*(exponent_next+exponent);
      PDF[i] += step.imag();
      
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda = lambda_next;
      tau = tau_next;
      exponent = exponent_next;
      if(j == 0){
        first_step = step;
      }
      j++;
    }while(abs(step/first_step) > 1.0e-5 || j < int(6.0/sigma_frac));
    
    PDF[i] /= constants::pi;
    
    if(i%1 == 0){
      cout << PDF[i] << setw(20);
      cout  << '\n';
    }
    
  }
  
  cout << "Done.\n";
  
}




extern "C" void return_PDF_2D(double* delta_values, double* PDF, int n_delta, double delta_min, double delta_max, double z, double R_in_comoving_Mpc, double L_in_comoving_Mpc, double f_NL, double var_NL_rescale, int LOS_modus, int index_of_universe){
  
  cout << "Computing phi_data:\n";
  cout.flush();
  
  vector<vector<double> > phi_data;
  
  phi_data = global_universes.matter_contents[index_of_universe]->compute_phi_of_lambda_2D(z, R_in_comoving_Mpc/constants::c_over_e5, L_in_comoving_Mpc/constants::c_over_e5, f_NL, var_NL_rescale, LOS_modus);
  
  /*
   * Determine critical point, where phi(lambda) splits into two branches on the complex plane. 
   * 
   */
  int n_lambda = 0;
  double lambda_c = phi_data[2][0];
  double delta_c = phi_data[1][0];
  double tau_c = phi_data[0][0]/sqrt(phi_data[6][0]);
  for(int i = 1; i < phi_data[4].size(); i++){
    if(phi_data[4][i-1] < phi_data[4][i]){
      n_lambda = i+1;
      lambda_c = phi_data[4][i];
      delta_c = phi_data[1][i];
      tau_c = phi_data[0][i]/sqrt(phi_data[6][i]);
    }
    else{
      i = 2*phi_data[4].size();
    }
  }
  
  /*
   * Extract phi_data up to the critical point.
   * 
   */
  vector<double> delta_L_values(n_lambda, 0.0);
  vector<double> delta_NL_values(n_lambda, 0.0);
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> phi_values(n_lambda, 0.0);
  vector<double> lambda_values_Gauss(n_lambda, 0.0);
  vector<double> phi_values_Gauss(n_lambda, 0.0);
  vector<double> variances(n_lambda, 0.0);
  vector<double> skewnesses(n_lambda, 0.0);
  vector<double> R_L_values(n_lambda, 0.0);
  vector<double> tau_values(n_lambda, 0.0);
  
  for(int i = 0; i < tau_values.size(); i++){
    delta_L_values[i] = phi_data[0][i];
    delta_NL_values[i] = phi_data[1][i];
    lambda_values[i] = phi_data[2][i];
    phi_values[i] = phi_data[3][i];
    lambda_values_Gauss[i] = phi_data[4][i];
    phi_values_Gauss[i] = phi_data[5][i];
    variances[i] = phi_data[6][i];
    skewnesses[i] = phi_data[7][i];
    R_L_values[i] = phi_data[8][i];
    tau_values[i] = delta_L_values[i]/sqrt(variances[i]);
  }
  
  
  /*
   * Extract phi_data with equal number and range of points left and right of tau = 0 (better for polynomial fit).
   * 
   */
  
  double tau_max = 0.9*tau_c;
  double tau_min = 0.9*tau_values[0];
  
  //double R_max = interpolate_neville_aitken(tau_max, &tau_values, &R_L_values, constants::order_of_interpolation);
  //if(R_max > exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5){
  //  tau_max = min(tau_max, interpolate_neville_aitken(exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  //}
  //double R_min = interpolate_neville_aitken(tau_min, &tau_values, &R_L_values, constants::order_of_interpolation);
  //if(R_min < exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5){
  //  tau_min = max(tau_min, interpolate_neville_aitken(exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  //}
  
  tau_c = tau_max;
  
  delta_c = interpolate_neville_aitken(tau_max, &tau_values, &delta_NL_values, constants::order_of_interpolation);
  lambda_c = interpolate_neville_aitken(tau_max, &tau_values, &lambda_values, constants::order_of_interpolation);
  
  int n_tau = 4*constants::generating_function_coeff_order + 1; // has to be odd number in order to include tau=0 exactly.
  
  vector<double> tau_for_fit(n_tau,0.0);
  vector<double> lambda_for_fit(n_tau,0.0);
  vector<double> phi_for_fit(n_tau,0.0);
  vector<double> phi_prime_for_fit(n_tau,0.0);
  
  
  double dt = -tau_min/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i] = tau_min+double(i)*dt;
  }
  tau_for_fit[n_tau/2] = 0.0;
  dt = tau_max/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i+n_tau/2+1] = double(i+1)*dt;
  }
  
  for(int i = 0; i < n_tau; i++){
    lambda_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &lambda_values, constants::order_of_interpolation);
    phi_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &phi_values, constants::order_of_interpolation);
    phi_prime_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &delta_NL_values, constants::order_of_interpolation);
  }
  
  cout << "Done.\n";
  
  /*
   * Express functions as polynomials in tau.
   * 
   */
  cout << "Computing tau coefficients:\n";
  
  int n_coeff = constants::generating_function_coeff_order;
  
  vector<double> coefficients_lambda_of_tau = return_coefficients(&tau_for_fit, &lambda_for_fit, n_coeff);
  vector<double> coefficients_lambda_of_tau_prime(coefficients_lambda_of_tau.size(), 0.0);
  for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++) coefficients_lambda_of_tau_prime[i] = coefficients_lambda_of_tau[i+1]*double(i+1);
  
  vector<double> coefficients_phi_of_tau = return_coefficients(&tau_for_fit, &phi_for_fit, n_coeff);
  vector<double> coefficients_phi_of_tau_prime = return_coefficients(&tau_for_fit, &phi_prime_for_fit, n_coeff);
  
  cout << "Done.\n";
  
  /*
   * Perform the inverse Laplace transform of phi(lambda) to compute p(delta).
   * 
   */
  double ddelta = (delta_max-delta_min)/double(n_delta-1);
  double delta;
  double tau_0, lambda_0;
  double dr;
  
  complex<double> lambda;
  complex<double> lambda_next;
  complex<double> tau, tau_next;
  complex<double> phi_prime, phi_prime_next;
  complex<double> exponent, exponent_next;
  complex<double> dlambda;
  complex<double> dlambda_dtau;
  complex<double> step, first_step;
  
  cout << "Computing PDF:\n";
  for(int i = 0; i < n_delta; i++){
    delta = delta_min + double(i)*ddelta;
    
    delta_values[i] = delta;
    PDF[i] = 0.0;
      
    if(delta < delta_c){
      tau_0 = interpolate_Newton(delta, &delta_NL_values, &tau_values, constants::order_of_interpolation);
      lambda_0 = interpolate_Newton(delta, &delta_NL_values, &lambda_values, constants::order_of_interpolation);
    }
    else{
      tau_0 = tau_c;
      lambda_0 = lambda_c;
    }
    lambda = complex<double>(lambda_0, 0.0);
    tau = complex<double>(tau_0, 0.0);
    exponent = exp(-lambda*delta + return_polnomial_value(tau, &coefficients_phi_of_tau));
    
    // sigma_r^2 \approx 1/phi''(lambda_0)
    double sigma_frac = 0.001;
    dr = sigma_frac/sqrt(interpolate_neville_aitken_derivative(lambda_0, &lambda_values, &delta_NL_values, constants::order_of_interpolation));
    dlambda = complex<double>(0.0, dr);
    int j = 0;
    
    if(i%1 == 0){
      cout << delta << setw(20);
      cout.flush();
    }
    do{
      lambda_next = lambda + 0.5*dlambda;
      tau_next = Newtons_method_complex(tau, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda_next = lambda + dlambda;
      tau_next = Newtons_method_complex(tau_next, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      exponent_next = exp(-lambda_next*delta + return_polnomial_value(tau_next, &coefficients_phi_of_tau));
      
      step = 0.5*dlambda*(exponent_next+exponent);
      PDF[i] += step.imag();
      
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda = lambda_next;
      tau = tau_next;
      exponent = exponent_next;
      if(j == 0){
        first_step = step;
      }
      j++;
    }while(abs(step/first_step) > 1.0e-5 || j < int(6.0/sigma_frac));
    
    PDF[i] /= constants::pi;
    
    if(i%1 == 0){
      cout << PDF[i] << setw(20);
      cout  << '\n';
    }
    
  }
  
  cout << "Done.\n";
  
}



/**************** CODE FOR EdS APPROXIMATION OF SPHERICAL COLLAPSE **********************/




extern "C" void return_CGF_EdS(double* delta_L, double* delta_NL, double* lambda, double* phi, double* lambda_Gauss, double* phi_Gauss, double* variances, double* skewnesses, double* R_L, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  vector<vector<double> > phi_data = global_universes.matter_contents[index_of_universe]->compute_phi_of_lambda_3D_EdS(z, R_in_comoving_Mpc/constants::c_over_e5, f_NL, var_NL_rescale);

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

extern "C" void return_PDF_EdS(double* delta_values, double* PDF, int n_delta, double delta_min, double delta_max, double z, double R_in_comoving_Mpc, double f_NL, double var_NL_rescale, int index_of_universe){
  
  cout << "Computing phi_data:\n";
  cout.flush();
  
  vector<vector<double> > phi_data;
  
  phi_data = global_universes.matter_contents[index_of_universe]->compute_phi_of_lambda_3D_EdS(z, R_in_comoving_Mpc/constants::c_over_e5, f_NL, var_NL_rescale);
  
  /*
   * Determine critical point, where phi(lambda) splits into two branches on the complex plane. 
   * 
   */
  int n_lambda = 0;
  double lambda_c = phi_data[2][0];
  double delta_c = phi_data[1][0];
  double tau_c = phi_data[0][0]/sqrt(phi_data[6][0]);
  for(int i = 1; i < phi_data[4].size(); i++){
    if(phi_data[4][i-1] < phi_data[4][i]){
      n_lambda = i+1;
      lambda_c = phi_data[4][i];
      delta_c = phi_data[1][i];
      tau_c = phi_data[0][i]/sqrt(phi_data[6][i]);
    }
    else{
      i = 2*phi_data[4].size();
    }
  }
  
  /*
   * Extract phi_data up to the critical point.
   * 
   */
  vector<double> delta_L_values(n_lambda, 0.0);
  vector<double> delta_NL_values(n_lambda, 0.0);
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> phi_values(n_lambda, 0.0);
  vector<double> lambda_values_Gauss(n_lambda, 0.0);
  vector<double> phi_values_Gauss(n_lambda, 0.0);
  vector<double> variances(n_lambda, 0.0);
  vector<double> skewnesses(n_lambda, 0.0);
  vector<double> R_L_values(n_lambda, 0.0);
  vector<double> tau_values(n_lambda, 0.0);
  
  for(int i = 0; i < tau_values.size(); i++){
    delta_L_values[i] = phi_data[0][i];
    delta_NL_values[i] = phi_data[1][i];
    lambda_values[i] = phi_data[2][i];
    phi_values[i] = phi_data[3][i];
    lambda_values_Gauss[i] = phi_data[4][i];
    phi_values_Gauss[i] = phi_data[5][i];
    variances[i] = phi_data[6][i];
    skewnesses[i] = phi_data[7][i];
    R_L_values[i] = phi_data[8][i];
    tau_values[i] = delta_L_values[i]/sqrt(variances[i]);
  }
  
  
  /*
   * Extract phi_data with equal number and range of points left and right of tau = 0 (better for polynomial fit).
   * 
   */
  
  double tau_max = 0.9*tau_c;
  double tau_min = 0.9*tau_values[0];
  
  double R_max = interpolate_neville_aitken(tau_max, &tau_values, &R_L_values, constants::order_of_interpolation);
  if(R_max > exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5){
    tau_max = min(tau_max, interpolate_neville_aitken(exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses.size()-1])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  }
  double R_min = interpolate_neville_aitken(tau_min, &tau_values, &R_L_values, constants::order_of_interpolation);
  if(R_min < exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5){
    tau_min = max(tau_min, interpolate_neville_aitken(exp(global_universes.matter_contents[index_of_universe]->log_top_hat_radii_for_skewnesses[0])*constants::c_over_e5, &R_L_values, &tau_values, constants::order_of_interpolation));
  }
  
  tau_c = tau_max;
  
  delta_c = interpolate_neville_aitken(tau_max, &tau_values, &delta_NL_values, constants::order_of_interpolation);
  lambda_c = interpolate_neville_aitken(tau_max, &tau_values, &lambda_values, constants::order_of_interpolation);
  
  int n_tau = 4*constants::generating_function_coeff_order + 1; // has to be odd number in order to include tau=0 exactly.
  
  vector<double> tau_for_fit(n_tau,0.0);
  vector<double> lambda_for_fit(n_tau,0.0);
  vector<double> phi_for_fit(n_tau,0.0);
  vector<double> phi_prime_for_fit(n_tau,0.0);
  
  
  double dt = -tau_min/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i] = tau_min+double(i)*dt;
  }
  tau_for_fit[n_tau/2] = 0.0;
  dt = tau_max/double(n_tau/2);
  for(int i = 0; i < n_tau/2; i++){
    tau_for_fit[i+n_tau/2+1] = double(i+1)*dt;
  }
  
  for(int i = 0; i < n_tau; i++){
    lambda_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &lambda_values, constants::order_of_interpolation);
    phi_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &phi_values, constants::order_of_interpolation);
    phi_prime_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &delta_NL_values, constants::order_of_interpolation);
  }
  
  cout << "Done.\n";
  
  /*
   * Express functions as polynomials in tau.
   * 
   */
  cout << "Computing tau coefficients:\n";
  
  int n_coeff = constants::generating_function_coeff_order;
  
  vector<double> coefficients_lambda_of_tau = return_coefficients(&tau_for_fit, &lambda_for_fit, n_coeff);
  vector<double> coefficients_lambda_of_tau_prime(coefficients_lambda_of_tau.size(), 0.0);
  for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++) coefficients_lambda_of_tau_prime[i] = coefficients_lambda_of_tau[i+1]*double(i+1);
  
  vector<double> coefficients_phi_of_tau = return_coefficients(&tau_for_fit, &phi_for_fit, n_coeff);
  vector<double> coefficients_phi_of_tau_prime = return_coefficients(&tau_for_fit, &phi_prime_for_fit, n_coeff);
  
  cout << "Done.\n";
  
  /*
   * Perform the inverse Laplace transform of phi(lambda) to compute p(delta).
   * 
   */
  double ddelta = (delta_max-delta_min)/double(n_delta-1);
  double delta;
  double tau_0, lambda_0;
  double dr;
  
  complex<double> lambda;
  complex<double> lambda_next;
  complex<double> tau, tau_next;
  complex<double> phi_prime, phi_prime_next;
  complex<double> exponent, exponent_next;
  complex<double> dlambda;
  complex<double> dlambda_dtau;
  complex<double> step, first_step;
  
  cout << "Computing PDF:\n";
  for(int i = 0; i < n_delta; i++){
    delta = delta_min + double(i)*ddelta;
    
    delta_values[i] = delta;
    PDF[i] = 0.0;
      
    if(delta < delta_c){
      tau_0 = interpolate_Newton(delta, &delta_NL_values, &tau_values, constants::order_of_interpolation);
      lambda_0 = interpolate_Newton(delta, &delta_NL_values, &lambda_values, constants::order_of_interpolation);
    }
    else{
      tau_0 = tau_c;
      lambda_0 = lambda_c;
    }
    lambda = complex<double>(lambda_0, 0.0);
    tau = complex<double>(tau_0, 0.0);
    exponent = exp(-lambda*delta + return_polnomial_value(tau, &coefficients_phi_of_tau));
    
    // sigma_r^2 \approx 1/phi''(lambda_0)
    double sigma_frac = 0.001;
    dr = sigma_frac/sqrt(interpolate_neville_aitken_derivative(lambda_0, &lambda_values, &delta_NL_values, constants::order_of_interpolation));
    dlambda = complex<double>(0.0, dr);
    int j = 0;
    
    if(i%1 == 0){
      cout << delta << setw(20);
      cout.flush();
    }
    do{
      lambda_next = lambda + 0.5*dlambda;
      tau_next = Newtons_method_complex(tau, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda_next = lambda + dlambda;
      tau_next = Newtons_method_complex(tau_next, lambda_next,  &coefficients_lambda_of_tau, &coefficients_lambda_of_tau_prime);
      phi_prime_next = return_polnomial_value(tau_next, &coefficients_phi_of_tau_prime);
      exponent_next = exp(-lambda_next*delta + return_polnomial_value(tau_next, &coefficients_phi_of_tau));
      
      step = 0.5*dlambda*(exponent_next+exponent);
      PDF[i] += step.imag();
      
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda = lambda_next;
      tau = tau_next;
      exponent = exponent_next;
      if(j == 0){
        first_step = step;
      }
      j++;
    }while(abs(step/first_step) > 1.0e-5 || j < int(6.0/sigma_frac));
    
    PDF[i] /= constants::pi;
    
    if(i%1 == 0){
      cout << PDF[i] << setw(20);
      cout  << '\n';
    }
    
  }
  
  cout << "Done.\n";
  
}



/************* END OF CODE FOR EdS APPROXIMATION OF SPHERICAL COLLAPSE ******************/



extern "C" int return_Nlambda(int index_of_universe){
  return global_universes.matter_contents[index_of_universe]->return_N_of_lambda();
}


extern "C" int return_Nlambda_EdS(int index_of_universe){
  vector<double> delta_L_values;
  vector<double> delta_NL_values;
  vector<double> delta_NL_prime_values;
  
  return_EdS_spherical_collapse(1.0, &delta_L_values, &delta_NL_values, &delta_NL_prime_values);
  
  return delta_L_values.size();
}


extern "C" double return_var_NL(double z, double R_in_Mpc_over_h, int index_of_universe){
  return global_universes.matter_contents[index_of_universe]->return_non_linear_variance(z, R_in_Mpc_over_h);
}


extern "C" double return_var_L(double z, double R_in_Mpc_over_h, int index_of_universe){
  return global_universes.matter_contents[index_of_universe]->return_linear_variance(z, R_in_Mpc_over_h);
}

extern "C" void return_var_NL_2D(double *var_deltaLOS, double *var_GaussLOS, double *var_tophatLOS, double z, double R_in_Mpc_over_h, double L_in_Mpc_over_h, int index_of_universe){
  global_universes.matter_contents[index_of_universe]->return_2D_non_linear_variance(var_deltaLOS, var_GaussLOS, var_tophatLOS, z, R_in_Mpc_over_h, L_in_Mpc_over_h);
}


