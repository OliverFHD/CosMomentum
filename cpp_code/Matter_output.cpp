
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

