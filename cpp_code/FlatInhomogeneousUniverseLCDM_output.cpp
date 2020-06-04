
/**********************************************
 **********************************************
 **_________ 5. OUTPUT AND CHECKS ___________**
 **********************************************
 ********************************************** 
 *                                            *
 * ..... 5.1 print_Newtonian_growth_factor    *
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

void FlatInhomogeneousUniverseLCDM::print_Newtonian_growth_factor(string file_name){
  
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
  
  for(int i = 0; i < this->return_number_of_time_steps(); i++){
    e = this->eta[i];
    scale = this->a[i];
    hubble = this->H[i];
    z = 1.0/scale - 1.0;
    D = this->Newtonian_growth_factor_of_delta[i];
    D_prime = this->Newtonian_growth_factor_of_delta_prime[i];
    output << setw(20) << e << setw(20) << scale << setw(20) << hubble << setw(20) << z << setw(20) << D << setw(20) << D_prime << setw(20) << D/scale << setw(20) << D_prime/scale - D/scale*hubble << '\n';
  }
  
}

/*******************************************************************************************************************************************************
 * 5.3 transfer_function_at 
 * Description:
 *
 * Arguments:
 * 
 * 
*******************************************************************************************************************************************************/

double FlatInhomogeneousUniverseLCDM::transfer_function_at(double k){
  return interpolate_neville_aitken(log(k), &this->log_wave_numbers, &this->transfer_function, constants::order_of_interpolation);
}


void FlatInhomogeneousUniverseLCDM::return_delta_NL_of_delta_L_and_dF_ddelta_3D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values){
  (*delta_L_values) = this->delta_values_for_spherical_collapse;
  int N_delta = this->delta_values_for_spherical_collapse.size();
  (*delta_NL_values) = vector<double>(N_delta);
  (*delta_NL_prime_values) = vector<double>(N_delta);
  for(int i = 0; i < N_delta; i++){
    (*delta_NL_values)[i] = interpolate_neville_aitken(eta, &this->eta, &this->spherical_collapse_evolution_of_delta[i], constants::order_of_interpolation);
    (*delta_NL_prime_values)[i] = interpolate_neville_aitken(eta, &this->eta, &this->spherical_collapse_evolution_of_delta_ddelta[i], constants::order_of_interpolation);
  }
}


void FlatInhomogeneousUniverseLCDM::return_delta_NL_of_delta_L_and_dF_ddelta_2D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values){
  (*delta_L_values) = this->delta_values_for_cylindrical_collapse;
  int N_delta = this->delta_values_for_cylindrical_collapse.size();
  (*delta_NL_values) = vector<double>(N_delta);
  (*delta_NL_prime_values) = vector<double>(N_delta);
  for(int i = 0; i < N_delta; i++){
    (*delta_NL_values)[i] = interpolate_neville_aitken(eta, &this->eta, &this->cylindrical_collapse_evolution_of_delta[i], constants::order_of_interpolation);
    (*delta_NL_prime_values)[i] = interpolate_neville_aitken(eta, &this->eta, &this->cylindrical_collapse_evolution_of_delta_ddelta[i], constants::order_of_interpolation);
  }
}


double FlatInhomogeneousUniverseLCDM::return_D_of_z(double z){
  
  double e = this->eta_at_a(1.0/(1.0+z));
  return interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  
}

double FlatInhomogeneousUniverseLCDM::return_D_of_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  
}

vector<vector<double> > FlatInhomogeneousUniverseLCDM::return_linear_growth_history(int conformal_time_steps){
  
  int n_column = 4;
  vector<vector<double> > linear_growth_history(n_column, vector<double>(conformal_time_steps, 0.0));
  
  double e_max = this->return_eta_final();
  double e_min = this->return_eta_initial();
  double de = (e_max - e_min)/double(conformal_time_steps-1);
  double e, a;
  
  for(int i = 0; i < conformal_time_steps; i++){
    e = e_min + double(i)*de;
    a = this->a_at_eta(e);
    linear_growth_history[0][i] = e;
    linear_growth_history[1][i] = a;
    linear_growth_history[2][i] = this->return_D_of_eta(e);
    linear_growth_history[3][i] = this->return_D_prime_of_eta(e);
  }
  
  return linear_growth_history;
  
}

double FlatInhomogeneousUniverseLCDM::return_D_prime_of_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta_prime, constants::order_of_interpolation);
  
}

vector<vector<double> > FlatInhomogeneousUniverseLCDM::return_power_spectra(double e, double R){
  
  double D = this->return_D_of_eta(e);
  double k;
  
  this->current_P_L = this->P_L(e);
  this->current_P_NL = this->P_NL(e);
  vector<vector<double> > power_spectra(3, vector<double>(constants::number_of_k,0.0));
  
  for(int i = 0; i < constants::number_of_k; i++){
    k = this->wave_numbers[i];
    power_spectra[0][i] = k/c_over_e5;
    power_spectra[1][i] = this->current_P_L[i]*pow(c_over_e5,3.0);
    power_spectra[2][i] = this->current_P_NL[i]*pow(c_over_e5,3.0);
  }
  
  return power_spectra;
  
}

void FlatInhomogeneousUniverseLCDM::return_2nd_moment_and_derivative(double R, double *variance, double *dvariance_dR){
  
  (*variance) = interpolate_neville_aitken(log(R), &this->log_top_hat_radii, &this->top_hat_sphere_variances, constants::order_of_interpolation);
  (*dvariance_dR) = interpolate_neville_aitken(log(R), &this->log_top_hat_radii, &this->dtop_hat_sphere_variances_dR, constants::order_of_interpolation);
  
}

void FlatInhomogeneousUniverseLCDM::return_2nd_moment_and_derivative_2D(double R, double *variance, double *dvariance_dR){
  
  (*variance) = interpolate_neville_aitken(log(R), &this->log_top_hat_cylinder_radii, &this->top_hat_cylinder_variances, constants::order_of_interpolation);
  (*dvariance_dR) = interpolate_neville_aitken(log(R), &this->log_top_hat_cylinder_radii, &this->dtop_hat_cylinder_variances_dR, constants::order_of_interpolation);
  
}

double FlatInhomogeneousUniverseLCDM::return_non_linear_variance(double z, double R_in_Mpc_over_h){
  
  double eta = this->eta_at_a(1.0/(1.0+z));
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  
  this->current_P_NL = this->P_NL(eta);
  double var_NL_R = variance_of_matter_within_R_NL(R);
  
  return var_NL_R;
  
}

double FlatInhomogeneousUniverseLCDM::return_linear_variance(double z, double R_in_Mpc_over_h){
  
  double eta = this->eta_at_a(1.0/(1.0+z));
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  
  this->current_P_L = this->P_L(eta);
  double var_L_R = variance_of_matter_within_R(R);
  
  return var_L_R;
  
}


void FlatInhomogeneousUniverseLCDM::print_growth_history(string file_name){
  
  double Om_m = this->return_Omega_m();
  double Om_l = this->return_Omega_L();
  double Om_r = this->return_Omega_r();
  double scale, e;
  fstream growth_stream;
  
  remove(file_name.c_str());
  FILE * F = fopen(file_name.c_str(), "w");
  fclose(F);
  growth_stream.open(file_name.c_str());
  
  growth_stream << setprecision(10) << scientific;
  growth_stream << "#Cosmological Parameters: Om_m = " << Om_m << ", Om_l = " << Om_l << ", Om_r = " << Om_r << '\n';
  growth_stream << "#a_max = " << this->a_final << '\n';
  growth_stream << "#eta(t)" << setw(20) << "w(eta)" << setw(20) << "z(eta)" << setw(20) << "a(eta)" << setw(20) << "D(eta)" << setw(20) << "D_prime(eta)\n";
  for(int i = 0; i < this->return_number_of_time_steps(); i++){
    e = this->eta[i];
    scale = this->a[i];
    growth_stream << setw(20) << e;
    growth_stream << setw(20) << this->eta_at_a(1.0)-e;
    growth_stream << setw(20) << 1.0/scale - 1.0;
    growth_stream << setw(20) << scale;
    growth_stream << setw(20) << this->Newtonian_growth_factor_of_delta[i];
    growth_stream << setw(20) << this->Newtonian_growth_factor_of_delta_prime[i] << '\n';
  }
    
  growth_stream.close();
  
}

/*******************************************************************************************************************************************************
 * return_LOS_integrated_phi_of_lambda
 * Description:
 * - set phi(delta, eta) and lambda(delta, eta) on a 2D grid. This grid is used when computing the LOS-projected CGF in Limber approximation
 * Arguments:
 * - double theta: angular top-hat radius with which LOS-integrated field is smoothed
 * 
*******************************************************************************************************************************************************/

vector<vector<double> > FlatInhomogeneousUniverseLCDM::return_LOS_integrated_phi_of_lambda(double theta, double f_NL, vector<double> w_values, vector<double> kernel_values){
  
  int n_lambda = this->delta_values_for_cylindrical_collapse.size();
  int n_time = w_values.size() - 1;
  double eta;
  double eta_0 = this->eta_at_a(1.0);
  double w;
  double dw;
  double R;
  
  vector<double> w_values_bin_center(n_time, 0.0);
  for(int i = 0; i < n_time; i++){
    w_values_bin_center[i] = 0.5*(w_values[i+1]+w_values[i]);
  }
  vector<int> indeces_of_nonzero_kernel(0,0);
  for(int i = 0; i < n_time; i++){
    if(kernel_values[i] != 0.0){
      indeces_of_nonzero_kernel.push_back(i);
    }
  }
  
  int n_time_of_nonzero_kernel = indeces_of_nonzero_kernel.size();
  vector<vector<double> > y_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_prime_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
  vector<vector<double> > dummy_data;
  
  double y_min = 0.0;
  double y_max = 0.0;
  vector<double>::iterator y_min_iterator;
  vector<double>::iterator y_max_iterator;
  int time_index_y_max, lambda_index_y_max, lambda_index_y_min;
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time_of_nonzero_kernel; t++){
    int i = indeces_of_nonzero_kernel[t];
    w = w_values_bin_center[i];
    eta = eta_0-w;
    R = w*theta;
    cout.flush();
    dummy_data = compute_phi_tilde_of_lambda_2D(eta, R, f_NL);
    // ISSUE: this order of indeces in confusing.
    y_values[t] = dummy_data[2];
    phi_values[t] = dummy_data[3];
    phi_prime_values[t] = dummy_data[1];
    for(int l = 0; l < n_lambda; l++){
      y_values[t][l] /= kernel_values[i];
    }
    
    // Determining the boundaries of y over which we can perform the projection intergral.
    // Also, finding the maximal lambda to which primary branch of CGF(\lambda) extends.
    y_max_iterator = std::max_element(y_values[t].begin(), y_values[t].end());
    if(t == 0){
      y_max=*y_max_iterator;
      time_index_y_max = t;
    }
    else if(*y_max_iterator<y_max){
      y_max=*y_max_iterator;
      time_index_y_max = t;
    }
    
    y_min_iterator = std::min_element(y_values[t].begin(), y_values[t].end());
    if(t == 0)
      y_min=*y_min_iterator;
    else if(*y_min_iterator>y_min)
      y_min=*y_min_iterator;
    
    lambda_index_y_max = std::distance(y_values[t].begin(), y_max_iterator);
    cout << phi_prime_values[t][lambda_index_y_max] << '\n';
    y_values[t] = vector<double>(&y_values[t][0], &y_values[t][lambda_index_y_max]+1);
    phi_values[t] = vector<double>(&phi_values[t][0], &phi_values[t][lambda_index_y_max]+1);
    phi_prime_values[t] = vector<double>(&phi_prime_values[t][0], &phi_prime_values[t][lambda_index_y_max]+1);
    
  }
  
  lambda_index_y_min = 0;
  lambda_index_y_max = y_values[time_index_y_max].size()-1;
  while(y_values[time_index_y_max][lambda_index_y_min] < y_min){
    lambda_index_y_min++;
  }
  
  int n_lambda_1st_branch = 1 + lambda_index_y_max - lambda_index_y_min;
  vector<vector<double> > projected_phi_data(4, vector<double>(n_lambda_1st_branch, 0.0));
  projected_phi_data[0] = vector<double>(&y_values[time_index_y_max][lambda_index_y_min], &y_values[time_index_y_max][lambda_index_y_max]+1);
  
  cout << "projecting CGF\n";
  
  for(int y = 0; y < n_lambda_1st_branch; y++){
    for(int t = 0; t < n_time_of_nonzero_kernel; t++){
      int i = indeces_of_nonzero_kernel[t];
      dw = w_values[i+1]-w_values[i];
      projected_phi_data[1][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &phi_values[t], constants::order_of_interpolation)*dw;
      projected_phi_data[2][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &phi_prime_values[t], constants::order_of_interpolation)*kernel_values[i]*dw;
      projected_phi_data[3][y] += interpolate_neville_aitken_derivative(projected_phi_data[0][y], &y_values[t], &phi_prime_values[t], constants::order_of_interpolation)*pow(kernel_values[i], 2)*dw;
    }
  }
  
  cout << "delta_min = " << projected_phi_data[2][0] << '\n';
  cout << "delta_max = " << projected_phi_data[2][n_lambda_1st_branch-1] << '\n';
  
  return projected_phi_data;
  
}


/*******************************************************************************************************************************************************
 * return_LOS_integrated_variance
 * Description:
 * - set phi(delta, eta) and lambda(delta, eta) on a 2D grid. This grid is used when computing the LOS-projected CGF in Limber approximation
 * Arguments:
 * - double theta: angular top-hat radius with which LOS-integrated field is smoothed
 * 
*******************************************************************************************************************************************************/

double FlatInhomogeneousUniverseLCDM::return_LOS_integrated_variance(double theta, vector<double> w_values, vector<double> kernel_values){
  
  int n_time = w_values.size()-1;
  
  vector<int> indeces_of_nonzero_kernel(0,0);
  
  double a, eta, eta_0, w, dw, R;
  eta_0 = this->eta_at_a(1.0);
  
  for(int i = 0; i < n_time; i++){
    if(kernel_values[i] > 0.0){
      indeces_of_nonzero_kernel.push_back(i);
    }
  }
  
  int n_time_of_nonzero_kernel = indeces_of_nonzero_kernel.size();
  double var_projected = 0.0;
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time_of_nonzero_kernel; t++){
    int i = indeces_of_nonzero_kernel[t];
    w = 0.5*(w_values[i+1]+w_values[i]);
    dw = w_values[i+1]-w_values[i];
    eta = eta_0-w;
    R = w*theta;
    this->current_P_NL = this->P_NL(eta);
    var_projected += variance_of_matter_within_R_NL_2D(R)*pow(kernel_values[i], 2)*dw;
  }
  
  return var_projected;
  
}


/*******************************************************************************************************************************************************
 * return_LOS_integrated_skewness
 * Description:
 * - set phi(delta, eta) and lambda(delta, eta) on a 2D grid. This grid is used when computing the LOS-projected CGF in Limber approximation
 * Arguments:
 * - double theta: angular top-hat radius with which LOS-integrated field is smoothed
 * 
*******************************************************************************************************************************************************/

double FlatInhomogeneousUniverseLCDM::return_LOS_integrated_skewness(double theta, double f_NL, vector<double> w_values, vector<double> kernel_values){
  
  if(f_NL != 0.0){
    cerr << "CAREFUL: f_NL != 0 not yet implemented in return_LOS_integrated_skewness.\n";
  }
  
  int n_time = w_values.size()-1;
  double a, e, eta_0, w, dw, R;
  eta_0 = this->eta_at_a(1.0);
  
  vector<int> indeces_of_nonzero_kernel(0,0);
  for(int i = 0; i < n_time; i++){
    if(kernel_values[i] > 0.0){
      indeces_of_nonzero_kernel.push_back(i);
    }
  }
  
  int n_time_of_nonzero_kernel = indeces_of_nonzero_kernel.size();
  double skewness_projected = 0.0;
  
  double D_11;
  double D_22;
  double mu;
  double one_plus_mu;
  double vNL, vL, dlnvL_dlnR, S_3;
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time_of_nonzero_kernel; t++){
    int i = indeces_of_nonzero_kernel[t];
    w = 0.5*(w_values[i+1]+w_values[i]);
    dw = w_values[i+1]-w_values[i];
    e = eta_0-w;
    R = w*theta;
    this->current_P_L = this->P_L(e);
    this->current_P_NL = this->P_NL(e);
    vNL = variance_of_matter_within_R_NL_2D(R);
    vL = this->variance_of_matter_within_R_2D(R);
    dlnvL_dlnR = this->dvariance_of_matter_within_R_dR_2D(R)/vL*R;
    
    D_11 = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
    D_22 = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_second_order, constants::order_of_interpolation);
    mu = 1.0 - D_22/D_11/D_11;
    one_plus_mu = (1.0+mu);
    
    S_3 = 3.0*one_plus_mu + 1.5*dlnvL_dlnR;
    skewness_projected += S_3*vNL*vNL*pow(kernel_values[i], 3)*dw;
  }
  
  return skewness_projected;
  
}

double FlatInhomogeneousUniverseLCDM::return_3D_skewness(double z, double R_in_Mpc_over_h, double f_NL){
  
  if(f_NL != 0.0){
    cerr << "CAREFUL: f_NL != 0 not yet implemented in return_3D_skewness.\n";
  }
  cerr << "CAREFUL: return_3D_skewness still uses EdS S_3.\n";
  
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  double e = this->eta_at_a(1.0/(1.0+z));
  
  double D_11 = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
  double D_22 = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_second_order, constants::order_of_interpolation);
  double mu = 1.0 - D_22/D_11/D_11;
  double one_plus_mu = (1.0+mu);
  
  this->current_P_L = this->P_L(e);
  this->current_P_NL = this->P_NL(e);
  
  double vNL = variance_of_matter_within_R_NL(R);
  double vL = this->variance_of_matter_within_R(R);
  double dlnvL_dlnR = this->dvariance_of_matter_within_R_dR(R)/vL*R;
  
  double S_3 = 34.0/7.0 + dlnvL_dlnR;
  
  return S_3*vNL*vNL;
  
}






/*
 * FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF
 * 
 * Returns PDF of line-of-sight projected matter density contrast with projection kernel specified by the input arrays w_values and kernel_values. w_values should store the lower-redshift edge of each histogram bin and n_of_w_values should contain the redshift histrogram (the histogram doesn't need to be normalised, since normalisation is enforced later on in the code). First column of the returned array is \delta smoothed with a spherical top-hat of R_in_Mpc_over_h, while 2nd column is p(\delta).
 * 
 */



vector<vector<double> > FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF(vector<double> w_values, vector<double> kernel_values, double theta, double f_NL, double var_NL_rescale){
  
  cout << "Computing projected phi_data:\n";
  cout.flush();
    
  vector<vector<double> > phi_data = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, w_values, kernel_values);
  
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
  
  
  /*
   * ISSUE: one can enforce even more coefficients to their analytical value!
   */
  coefficients_phi_of_tau[0] = 0.0;
  coefficients_phi_of_tau[1] = 0.0;
  coefficients_phi_of_tau_prime[0] = 0.0;
  
  cout << "Done.\n";
  
  /*
   * Perform the inverse Laplace transform of phi(lambda) to compute p(delta).
   * 
   */
  
  int n_delta = constants::N_delta_values_for_PDFs;
  double var = this->return_LOS_integrated_variance(theta, w_values, kernel_values);
  double delta_min = interpolate_neville_aitken(tau_min, &tau_values, &delta_NL_values, constants::order_of_interpolation);
  double delta_max = max(delta_c, 6.0*sqrt(var));
  // --> ISSUE: choosing delta_max to be 6*sigma may not be 100% reasonable for very skewed PDFs
  //            Maybe choose it by some quantile in a log-normal PDF that approximates the real PDF?
  
  
  /*
   * ISSUE: sometimes at delta=delta_max the code outputs PDF[i] = nan! Why is this? OF preliminarily fixed this by replacing "/double(n_delta-1)" with "/double(n_delta)" (hence avoiding delta=delta_max).
   */
  double ddelta = (delta_max-delta_min)/double(n_delta);
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

/*
 * FlatInhomogeneousUniverseLCDM::compute_polynomial_coefficients_from_CGF
 * 
 * Extract polynomial coefficients approximating the cumulant generating function. The problem here is: just fitting a polynomial to phi(lambda) is usually highly nummerically unstable. Instead, phi(tau) and lambda(tau) should be fitted by polynomials. From those polynomial coefficients one can then compute the polynomial coefficients of phi(lambda) (cf. section 4.6 of https://arxiv.org/pdf/1912.06621.pdf ).
 * 
 * Note: the nth polynomial coefficient a_n is connect to the nth cumulant c_n via
 * 
 * a_n = c_n / n!  .
 * 
 */

void FlatInhomogeneousUniverseLCDM::compute_polynomial_coefficients_from_CGF(vector<double> lambda_values, vector<double> tau_values, vector<double> phi_values, vector<double> phi_prime_values, vector<double> *coeffs_phi_of_lambda, vector<double> *coeffs_phi_prime_of_lambda){
  
  int n_lambda = lambda_values.size();
  
  int tau_c = tau_values[0];
  for(int i = 0; i < n_lambda-1; i++){
    if(lambda_values[i] < lambda_values[i+1]){
      tau_c = tau_values[i];
    }
  }
  
  
  double tau_max = 0.9*tau_c;
  double tau_min = 0.9*tau_values[0];
  int n_tau = 4*constants::generating_function_coeff_order + 1; // has to be odd number in order to include tau=0 exactly.
  // ISSUE: this shouldn't be hardcoded!
  
  vector<double> tau_for_fit(n_tau,0.0);
  vector<double> lambda_for_fit(n_tau,0.0);
  vector<double> phi_for_fit(n_tau,0.0);
  vector<double> phi_prime_for_fit(n_tau,0.0);
  
  // ISSUE: this only works as long as C++ doesn't change its treatment of integer division.
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
    phi_prime_for_fit[i] = interpolate_neville_aitken(tau_for_fit[i], &tau_values, &phi_prime_values, constants::order_of_interpolation);
  }
  
  // polynomial coefficients approximating lambda(tau)
  vector<double> coeffs_lambda_of_tau = return_coefficients(&tau_for_fit, &lambda_for_fit, constants::generating_function_coeff_order);
  int N_coeffs = coeffs_lambda_of_tau.size();
  
  // Bell matrices for translating coefficients interms of tau to coefficients in terms of lambda
  vector<vector<double> > Bell_matrix(0, vector<double>(0, 0.0));
  vector<vector<double> > inverse_Bell_matrix(0, vector<double>(0, 0.0));
  return_Bell_matrix(&Bell_matrix, &inverse_Bell_matrix, &coeffs_lambda_of_tau);
  
  
  vector<double> coeffs_phi_prime_of_tau = return_coefficients(&tau_for_fit, &phi_prime_for_fit, constants::order_of_interpolation);
  (*coeffs_phi_of_lambda) = vector<double>(N_coeffs,0.0);
  (*coeffs_phi_prime_of_lambda) = vector<double>(N_coeffs,0.0);
  vector<double> fact = factoria(N_coeffs);
  
  for(int i = 0; i < N_coeffs; i++){
    for(int j = 0; j <= i; j++){
      (*coeffs_phi_prime_of_lambda)[i] += inverse_Bell_matrix[i][j]*coeffs_phi_prime_of_tau[j]*fact[j];
    }
    (*coeffs_phi_prime_of_lambda)[i] /= fact[i];
  }
  
  
  for(int i = N_coeffs-1; i > 0; i--){
    (*coeffs_phi_of_lambda)[i] = (*coeffs_phi_prime_of_lambda)[i-1]/double(i);
  }
  
}





/*******************************************************************************************************************************************************
 * return_LOS_integrated_polynomial_coefficients
 * 
 * Returns polynomial coefficients of the cumulant generating function, projected along the line-of-sight kernel kernel_values.
 * NOTE: this function does not check whether the binning of w_values is too coarse for accurate integration.
 *       (The binning is instead set and controlled by ProjectedGalaxySample::set_n_of_w_data)
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::return_LOS_integrated_polynomial_coefficients(double theta, double f_NL, vector<double> w_values, vector<double> kernel_values, vector<double> *coeffs_phi_of_lambda, vector<double> *coeffs_phi_prime_of_lambda){
  
  int n_time = w_values.size()-1;
  double eta;
  double eta_0 = this->eta_at_a(1.0);
  double w;
  double dw;
  double R;
  
  vector<double> w_values_bin_center(n_time, 0.0);
  for(int i = 0; i < n_time; i++){
    w_values_bin_center[i] = 0.5*(w_values[i+1]+w_values[i]);
  }
  vector<int> indeces_of_nonzero_kernel(0,0);
  for(int i = 0; i < n_time; i++){
    if(kernel_values[i] != 0.0){
      indeces_of_nonzero_kernel.push_back(i);
    }
  }
  
  int n_time_of_nonzero_kernel = indeces_of_nonzero_kernel.size();
  vector<vector<double> > coefficients_phi_of_lambda(n_time_of_nonzero_kernel);
  vector<vector<double> > coefficients_phi_prime_of_lambda(n_time_of_nonzero_kernel);
  vector<vector<double> > CGF_data;
  vector<double> tau_values;
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time_of_nonzero_kernel; t++){
    int i = indeces_of_nonzero_kernel[t];
    w = w_values_bin_center[i];
    eta = eta_0-w;
    R = w*theta;
    
    CGF_data = this->compute_phi_tilde_of_lambda_2D(eta, R, f_NL);
    tau_values = CGF_data[0];
    for(int d = 0; d < tau_values.size(); d++){
      tau_values[d] /= sqrt(CGF_data[4][d]);
    }
    this->compute_polynomial_coefficients_from_CGF(CGF_data[2], tau_values, CGF_data[3], CGF_data[1], &coefficients_phi_of_lambda[t], &coefficients_phi_prime_of_lambda[t]);
    coefficients_phi_prime_of_lambda[t][0] = 0.0;
    // ISSUE: you can do this because current non-linear P(k) has been computed in "compute_phi_tilde_of_lambda_2D".
    // But this kind of flow seems very prone to error!
    coefficients_phi_prime_of_lambda[t][1] = variance_of_matter_within_R_NL_2D(R);
    coefficients_phi_of_lambda[t][0] = 0.0;
    coefficients_phi_of_lambda[t][1] = 0.0;
    coefficients_phi_of_lambda[t][2] = 0.5*coefficients_phi_prime_of_lambda[t][1];
  }
  
  int N_coeff = coefficients_phi_of_lambda[0].size();
  (*coeffs_phi_of_lambda) = vector<double>(N_coeff,0.0);
  (*coeffs_phi_prime_of_lambda) = vector<double>(N_coeff,0.0);
  
  for(int c = 0; c < N_coeff; c++){
    for(int t = 0; t < n_time_of_nonzero_kernel; t++){
      int i = indeces_of_nonzero_kernel[t];
      dw = (w_values[i+1]-w_values[i]);
      (*coeffs_phi_of_lambda)[c] += dw*coefficients_phi_of_lambda[t][c]*pow(kernel_values[i], c);
      (*coeffs_phi_prime_of_lambda)[c] += dw*coefficients_phi_prime_of_lambda[t][c]*pow(kernel_values[i], c+1);
    }
  }
  
}



/*******************************************************************************************************************************************************
 * return_LOS_integrated_polynomial_coefficients_incl_CMB_kappa
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::return_LOS_integrated_polynomial_coefficients_incl_CMB_kappa(double theta, double f_NL, vector<double> w_values, vector<double> kernel_values, vector<vector<double> > *coeffs_phi, vector<vector<double> > *coeffs_dphi_dlambda_delta, vector<vector<double> > *coeffs_dphi_dlambda_kappa){
  
  int n_lambda = this->delta_values_for_cylindrical_collapse.size();
  int n_time = w_values.size()-1;
  double z, a, eta, eta_0, w, dw, w_last_scattering, R;
  eta_0 = this->eta_at_a(1.0);
  w_last_scattering = eta_0-this->eta_at_a(1.0/(1.0+constants::z_last_scattering));
  
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
  double w_min = w_values_extended_to_last_scattering[0];
  double w_max = w_values_extended_to_last_scattering[n_time];
  double norm = 0.0;
  for(int i = 0; i < n_time; i++){
    w = w_values_bin_center[i];
    a = this->a_at_eta(eta_0-w);
    CMB_lensing_kernel[i] = 1.5*this->return_Omega_m()*w*(w_last_scattering-w)/w_last_scattering/a;
  }
  
  vector<vector<double> > coefficients_phi_of_lambda(n_time);
  vector<vector<double> > coefficients_phi_prime_of_lambda(n_time);
  vector<vector<double> > CGF_data;
  vector<double> tau_values;
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time; t++){
    w = w_values_bin_center[t];
    eta = eta_0-w;
    R = w*theta;
    
    CGF_data = this->compute_phi_tilde_of_lambda_2D(eta, R, f_NL);
    tau_values = CGF_data[0];
    for(int d = 0; d < tau_values.size(); d++){
      tau_values[d] /= sqrt(CGF_data[4][d]);
    }
    this->compute_polynomial_coefficients_from_CGF(CGF_data[2], tau_values, CGF_data[3], CGF_data[1], &coefficients_phi_of_lambda[t], &coefficients_phi_prime_of_lambda[t]);
    coefficients_phi_prime_of_lambda[t][0] = 0.0;
    // ISSUE: you can do this because current non-linear P(k) has been computed in "compute_phi_tilde_of_lambda_2D".
    // But this kind of flow seems very prone to error!
    coefficients_phi_prime_of_lambda[t][1] = variance_of_matter_within_R_NL_2D(R);
    coefficients_phi_of_lambda[t][0] = 0.0;
    coefficients_phi_of_lambda[t][1] = 0.0;
    coefficients_phi_of_lambda[t][2] = 0.5*coefficients_phi_prime_of_lambda[t][1];
  }
  
  int N_coeff = coefficients_phi_of_lambda[0].size();
  (*coeffs_phi) = vector<vector<double> >(N_coeff,vector<double>(N_coeff, 0.0));
  (*coeffs_dphi_dlambda_delta) = vector<vector<double> >(N_coeff,vector<double>(N_coeff, 0.0));
  (*coeffs_dphi_dlambda_kappa) = vector<vector<double> >(N_coeff,vector<double>(N_coeff, 0.0));
  
  for(int c_delta = 0; c_delta < N_coeff; c_delta++){
    for(int c_kappa = 0; c_kappa < N_coeff; c_kappa++){
      for(int t = 0; t < n_time; t++){
        if(c_delta + c_kappa < N_coeff){
          dw = w_values_extended_to_last_scattering[t+1]-w_values_extended_to_last_scattering[t];
          (*coeffs_phi)[c_delta][c_kappa] += dw*coefficients_phi_of_lambda[t][c_delta+c_kappa]*pow(kernel_values[t], c_delta)*pow(CMB_lensing_kernel[t], c_kappa);
          (*coeffs_dphi_dlambda_delta)[c_delta][c_kappa] += dw*coefficients_phi_prime_of_lambda[t][c_delta+c_kappa]*pow(kernel_values[t], c_delta+1)*pow(CMB_lensing_kernel[t], c_kappa);
          (*coeffs_dphi_dlambda_kappa)[c_delta][c_kappa] += dw*coefficients_phi_prime_of_lambda[t][c_delta+c_kappa]*pow(kernel_values[t], c_delta)*pow(CMB_lensing_kernel[t], c_kappa+1);
        }
      }
    }
  }
  
}






