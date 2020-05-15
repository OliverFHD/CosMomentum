
/**********************************************
 **********************************************
 **_________ 5. OUTPUT AND CHECKS ___________**
 **********************************************
 ********************************************** 
 *                                            *
 * ..... 5.1 print_Newtonian_growth_factor    *
 * ..... 5.2 return_wave_numbers              *
 * ..... 5.3 transfer_function_at             *
 *                                            *
 **********************************************
 **********************************************/



/*******************************************************************************************************************************************************
 * 5.1 print_Newtonian_growth_factor
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

void Matter::print_Newtonian_growth_factor(string file_name){
  
  double e;
  double scale;
  double hubble;
  double z;
  double D;
  double D_prime;
  double d_eta, dD;
  
  fstream output;
  
  remove(file_name.c_str());
  FILE* F = fopen(file_name.c_str(),"w");
  fclose(F);
  output.open(file_name.c_str());
  
  
  output << setw(20) << "eta" << setw(20) << "a" << setw(20) << "H_conf" << setw(20) << "z" << setw(20) << "D(z)" << setw(20) << "D'(z)" << setw(20) << "D_phi(z)" << setw(20) << "D_phi'(z)" << '\n';
  
  for(int i = 1; i < this->number_of_entries_Newton; i++){
    e = this->eta_Newton[i];
    d_eta = e - this->eta_Newton[i-1];
    scale = this->universe->a_at_eta(e);
    hubble = this->universe->H_at_eta(e);
    z = 1.0/scale - 1.0;
    D = this->Newtonian_growth_factor_of_delta[i];
    D_prime = this->Newtonian_growth_factor_of_delta_prime[i];
    dD = D - this->Newtonian_growth_factor_of_delta[i-1];
    output << setw(20) << e << setw(20) << scale << setw(20) << hubble << setw(20) << z << setw(20) << D << setw(20) << D_prime << setw(20) << D/scale << setw(20) << D_prime/scale - D/scale*hubble << '\n';
  }
  
}

/*******************************************************************************************************************************************************
 * 5.2 return_wave_numbers 
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

vector<double> Matter::return_wave_numbers(){
  return this->wave_numbers;
}



/*******************************************************************************************************************************************************
 * 5.3 transfer_function_at 
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double Matter::transfer_function_at(double k){
  return interpolate_neville_aitken(log(k), &this->log_wave_numbers, &this->transfer_function, constants::order_of_interpolation);
}


void Matter::return_delta_NL_of_delta_L_and_dF_ddelta_3D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values){
  (*delta_L_values) = this->delta_values_for_spherical_collapse;
  (*delta_NL_values) = this->delta_values_for_spherical_collapse;
  (*delta_NL_prime_values) = this->delta_values_for_spherical_collapse;
  for(int i = 0; i < this->delta_values_for_spherical_collapse.size(); i++){
    (*delta_NL_values)[i] = interpolate_neville_aitken(eta, &this->eta_NL_for_spherical_collapse, &this->spherical_collapse_evolution_of_delta[i], constants::order_of_interpolation);
    (*delta_NL_prime_values)[i] = interpolate_neville_aitken(eta, &this->eta_NL_for_spherical_collapse, &this->spherical_collapse_evolution_of_delta_ddelta[i], constants::order_of_interpolation);
  }
}


void Matter::return_delta_NL_of_delta_L_and_dF_ddelta_2D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values){
  (*delta_L_values) = this->delta_values_for_cylindrical_collapse;
  (*delta_NL_values) = this->delta_values_for_cylindrical_collapse;
  (*delta_NL_prime_values) = this->delta_values_for_cylindrical_collapse;
  for(int i = 0; i < this->delta_values_for_cylindrical_collapse.size(); i++){
    (*delta_NL_values)[i] = interpolate_neville_aitken(eta, &this->eta_NL_for_cylindrical_collapse, &this->cylindrical_collapse_evolution_of_delta[i], constants::order_of_interpolation);
    (*delta_NL_prime_values)[i] = interpolate_neville_aitken(eta, &this->eta_NL_for_cylindrical_collapse, &this->cylindrical_collapse_evolution_of_delta_ddelta[i], constants::order_of_interpolation);
  }
}


double Matter::return_D_of_z(double z){
  
  double eta = this->universe->eta_at_a(1.0/(1.0+z));
  return interpolate_neville_aitken(eta, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  
}

double Matter::return_D_of_eta(double eta){
  
  return interpolate_neville_aitken(eta, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  
}

vector<vector<double> > Matter::return_linear_growth_history(int conformal_time_steps){
  
  int n_column = 4;
  vector<vector<double> > linear_growth_history(n_column, vector<double>(conformal_time_steps, 0.0));
  
  double e_max = this->eta_Newton[this->eta_Newton.size()-1];
  double e_min = this->eta_Newton[0];
  double de = (e_max - e_min)/double(conformal_time_steps-1);
  double e, a;
  
  for(int i = 0; i < conformal_time_steps; i++){
    e = e_min + double(i)*de;
    a = this->universe->a_at_eta(e);
    linear_growth_history[0][i] = e;
    linear_growth_history[1][i] = a;
    linear_growth_history[2][i] = this->return_D_of_eta(e);
    linear_growth_history[3][i] = this->return_D_prime_of_eta(e);
  }
  
  return linear_growth_history;
  
}

double Matter::return_D_prime_of_eta(double eta){
  
  return interpolate_neville_aitken(eta, &this->eta_Newton, &this->Newtonian_growth_factor_of_delta_prime, constants::order_of_interpolation);
  
}

vector<vector<double> > Matter::return_power_spectra(double eta, double R){
  
  double D = this->return_D_of_eta(eta);
  double k;
  
  this->current_P_L = this->P_L(eta);
  this->current_P_NL = this->P_NL(eta);
  vector<vector<double> > power_spectra(3, vector<double>(constants::number_of_k,0.0));
  
  for(int i = 0; i < constants::number_of_k; i++){
    k = this->wave_numbers[i];
    power_spectra[0][i] = k/c_over_e5;
    power_spectra[1][i] = this->current_P_L[i]*pow(c_over_e5,3.0);
    power_spectra[2][i] = this->current_P_NL[i]*pow(c_over_e5,3.0);
  }
  
  return power_spectra;
  
}

void Matter::return_2nd_moment_and_derivative(double R, double *variance, double *dvariance_dR){
  
  (*variance) = interpolate_neville_aitken(log(R), &this->log_top_hat_radii, &this->top_hat_sphere_variances, constants::order_of_interpolation);
  (*dvariance_dR) = interpolate_neville_aitken(log(R), &this->log_top_hat_radii, &this->dtop_hat_sphere_variances_dR, constants::order_of_interpolation);
  
}

void Matter::return_2nd_moment_and_derivative_2D(double R, double *variance, double *dvariance_dR){
  
  (*variance) = interpolate_neville_aitken(log(R), &this->log_top_hat_cylinder_radii, &this->top_hat_cylinder_variances, constants::order_of_interpolation);
  (*dvariance_dR) = interpolate_neville_aitken(log(R), &this->log_top_hat_cylinder_radii, &this->dtop_hat_cylinder_variances_dR, constants::order_of_interpolation);
  
}

double Matter::return_non_linear_variance(double z, double R_in_Mpc_over_h){
  
  double eta = this->universe->eta_at_a(1.0/(1.0+z));
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  
  this->current_P_NL = this->P_NL(eta);
  double var_NL_R = variance_of_matter_within_R_NL(R);
  
  return var_NL_R;
  
}

void Matter::return_2D_non_linear_variance(double *var_deltaLOS, double *var_GaussLOS, double *var_tophatLOS, double z, double R_in_Mpc_over_h, double L_in_Mpc_over_h){
  double eta = this->universe->eta_at_a(1.0/(1.0+z));
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  double L = L_in_Mpc_over_h/constants::c_over_e5;
  
  this->current_P_NL = this->P_NL(eta);
  
  (*var_deltaLOS) = variance_of_matter_within_R_NL_2D(R)/L;
  (*var_GaussLOS) = variance_of_matter_within_R_NL_2D_GaussianLOS(R, L);
  (*var_tophatLOS) = variance_of_matter_within_R_NL_2D(R, L);
}

double Matter::return_linear_variance(double z, double R_in_Mpc_over_h){
  
  double eta = this->universe->eta_at_a(1.0/(1.0+z));
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  
  this->current_P_L = this->P_L(eta);
  double var_L_R = variance_of_matter_within_R(R);
  
  return var_L_R;
  
}


void Matter::print_growth_history(string file_name){
  
  double Om_m = this->cosmology.Omega_m;
  double Om_l = this->cosmology.Omega_L;
  double Om_r = this->cosmology.Omega_r;
  double a, e;
  fstream growth_stream;
  
  remove(file_name.c_str());
  FILE * F = fopen(file_name.c_str(), "w");
  fclose(F);
  growth_stream.open(file_name.c_str());
  
  growth_stream << setprecision(10) << scientific;
  growth_stream << "#Cosmological Parameters: Om_m = " << Om_m << ", Om_l = " << Om_l << ", Om_r = " << Om_r << '\n';
  growth_stream << "#a_max = " << this->a_final << '\n';
  growth_stream << "#eta(t)" << setw(20) << "w(eta)" << setw(20) << "z(eta)" << setw(20) << "a(eta)" << setw(20) << "D(eta)" << setw(20) << "D_prime(eta)\n";
  for(int i = 0; i < this->eta_Newton.size(); i++){
    e = this->eta_Newton[i];
    a = this->universe->a_at_eta(e);
    growth_stream << setw(20) << e;
    growth_stream << setw(20) << this->universe->eta_at_a(1.0)-e;
    growth_stream << setw(20) << 1.0/a - 1.0;
    growth_stream << setw(20) << a;
    growth_stream << setw(20) << this->Newtonian_growth_factor_of_delta[i];
    growth_stream << setw(20) << this->Newtonian_growth_factor_of_delta_prime[i] << '\n';
  }
    
  growth_stream.close();
  
}






/*******************************************************************************************************************************************************
 * return_LOS_integrated_variance
 * Description:
 * - set phi(delta, eta) and lambda(delta, eta) on a 2D grid. This grid is used when computing the LOS-projected CGF in Limber approximation
 * Arguments:
 * - double theta: angular top-hat radius with which LOS-integrated field is smoothed
 * 
*******************************************************************************************************************************************************/

double Matter::return_LOS_integrated_variance(double theta, vector<double> z_values, vector<double> n_of_z_values){
  
  int n_time = z_values.size()-1;
  
  vector<double> w_values(n_time+1, 0.0);
  vector<double> dw_values(n_time, 0.0);
  vector<double> dz_values(n_time, 0.0);
  vector<double> n_of_w_values(n_time, 0.0);
  vector<int> indeces_of_nonzero_nofw(0,0);
  
  double a, eta, eta_0, w, R;
  eta_0 = this->universe->eta_at_a(1.0);
  
  cout << "normalising pofz data\n";
  
  for(int i = 0; i < n_time+1; i++){
    a = 1.0/(1.0+z_values[i]);
    w_values[i] = eta_0 - this->universe->eta_at_a(a);
  }
  
  double norm = 0.0;
  for(int i = 0; i < n_time; i++){
    dz_values[i] = z_values[i+1]-z_values[i];
    dw_values[i] = w_values[i+1]-w_values[i];
    n_of_w_values[i] = n_of_z_values[i]*dz_values[i]/dw_values[i];
    norm += n_of_z_values[i]*dz_values[i];
    if(n_of_w_values[i] > 0.0){
      indeces_of_nonzero_nofw.push_back(i);
    }
  }
  for(int i = 0; i < n_time; i++){
    n_of_w_values[i] /= norm;
  }
  // ISSUE --> if n_of_w_values is supposed to represent a lensing kernel, then it shouldn't be normalised.
  
  int n_time_of_nonzero_nofw = indeces_of_nonzero_nofw.size();
  double var_projected = 0.0;
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time_of_nonzero_nofw; t++){
    int i = indeces_of_nonzero_nofw[t];
    w = 0.5*(w_values[i+1]+w_values[i]);
    eta = eta_0-w;
    R = w*theta;
    this->current_P_NL = this->P_NL(eta);
    var_projected += variance_of_matter_within_R_NL_2D(R)*pow(n_of_w_values[i], 2)*dw_values[i];
  }
  
  return var_projected;
  
}








/*
 * Matter::compute_LOS_projected_PDF
 * 
 * Returns PDF of line-of-sight projected matter density contrast with projection kernel specified by the input arrays z_values and n_of_z_values. z_values should store the lower-redshift edge of each histogram bin and n_of_z_values should contain the redshift histrogram (the histogram doesn't need to be normalised, since normalisation is enforced later on in the code). First column of the returned array is \delta smoothed with a spherical top-hat of R_in_Mpc_over_h, while 2nd column is p(\delta).
 * 
 */



vector<vector<double> > Matter::compute_LOS_projected_PDF(vector<double> z_values, vector<double> n_of_z_values, double theta, double f_NL, double var_NL_rescale){
  
  cout << "Computing projected phi_data:\n";
  cout.flush();
    
  vector<vector<double> > phi_data = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, z_values, n_of_z_values);
  
  /*
   * Determine critical point, where phi(lambda) splits into two branches on the complex plane. 
   * 
   */
  int n_lambda = 0;
  double lambda_c = phi_data[0][0];
  double delta_c = phi_data[2][0];
  double tau_c = 0.0;
  for(int i = 1; i < phi_data[0].size(); i++){
    if(phi_data[0][i-1] < phi_data[0][i]){
      n_lambda = i+1;
      lambda_c = phi_data[0][i];
      delta_c = phi_data[2][i];
      if(phi_data[2][i]*phi_data[0][i] > phi_data[1][i])
        tau_c = sqrt(2.0*(phi_data[2][i]*phi_data[0][i] - phi_data[1][i]));
    }
    else{
      i = 2*phi_data[4].size();
    }
  }
  
  /*
   * Extract phi_data up to the critical point.
   * 
   */
  vector<double> delta_NL_values(n_lambda, 0.0);
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> phi_values(n_lambda, 0.0);
  vector<double> tau_values(n_lambda, 0.0);
  
  for(int i = 0; i < tau_values.size(); i++){
    lambda_values[i] = phi_data[0][i];
    phi_values[i] = phi_data[1][i];
    delta_NL_values[i] = phi_data[2][i];
    if(phi_data[2][i]*phi_data[0][i] > phi_data[1][i]){
      tau_values[i] = sqrt(2.0*(phi_data[2][i]*phi_data[0][i] - phi_data[1][i]));
      if(lambda_values[i] < 0.0) tau_values[i] *= -1.0;
    }
  }
  
  
  /*
   * Extract phi_data with equal number and range of points left and right of tau = 0 (better for polynomial fit).
   * 
   */
  
  double tau_max = 0.9*tau_c;
  double tau_min = 0.9*tau_values[0];
  
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
  
  int n_delta = constants::N_delta_values_for_PDFs;
  double var = this->return_LOS_integrated_variance(theta, z_values, n_of_z_values);
  double delta_min = interpolate_neville_aitken(tau_min, &tau_values, &delta_NL_values, constants::order_of_interpolation);
  double delta_max = max(delta_c, 6.0*sqrt(var));
  // --> ISSUE: choosing delta_max to be 6*sigma may not be 100% reasonable for very skewed PDFs
  //            Maybe choose it by some quantile in a log-normal PDF that approximates the real PDF?
  
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
  
  vector<vector<double> > PDF_data(2, vector<double>(n_delta));
  
  cout << "Computing 2D PDF:\n";
  for(int i = 0; i < n_delta; i++){
    delta = delta_min + double(i)*ddelta;
    
    PDF_data[0][i] = delta;
    PDF_data[1][i] = 0.0;
      
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
      PDF_data[1][i] += step.imag();
      
      dlambda = -dr*conj(phi_prime_next-delta)/abs(phi_prime_next-delta);
      lambda = lambda_next;
      tau = tau_next;
      exponent = exponent_next;
      if(j == 0){
        first_step = step;
      }
      j++;
    }while(abs(step/first_step) > 1.0e-5 || j < int(6.0/sigma_frac));
    
    PDF_data[1][i] /= constants::pi;
    
    cout << PDF_data[0][i] << "   ";
    cout << PDF_data[1][i] << "\n";
    
  }
  cout << "Done.\n";
  
  return PDF_data;
  
}
