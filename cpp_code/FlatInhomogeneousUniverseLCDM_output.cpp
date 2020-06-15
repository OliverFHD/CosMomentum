
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


void FlatInhomogeneousUniverseLCDM::return_delta_NL_of_delta_L_and_dF_ddelta_2D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values, vector<double> *delta_NL_prime_prime_values){
  (*delta_L_values) = this->delta_values_for_cylindrical_collapse;
  int N_delta = this->delta_values_for_cylindrical_collapse.size();
  (*delta_NL_values) = vector<double>(N_delta);
  (*delta_NL_prime_values) = vector<double>(N_delta);
  (*delta_NL_prime_prime_values) = vector<double>(N_delta);
  for(int i = 0; i < N_delta; i++){
    (*delta_NL_values)[i] = interpolate_neville_aitken(eta, &this->eta, &this->cylindrical_collapse_evolution_of_delta[i], constants::order_of_interpolation);
    (*delta_NL_prime_values)[i] = interpolate_neville_aitken(eta, &this->eta, &this->cylindrical_collapse_evolution_of_delta_ddelta[i], constants::order_of_interpolation);
    (*delta_NL_prime_prime_values)[i] = interpolate_neville_aitken(eta, &this->eta, &this->cylindrical_collapse_evolution_of_delta_ddelta2[i], constants::order_of_interpolation);
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

double FlatInhomogeneousUniverseLCDM::return_non_linear_variance(double z, double R_in_Mpc_over_h, double var_NL_rescale){
  
  double eta = this->eta_at_a(1.0/(1.0+z));
  double R = R_in_Mpc_over_h/constants::c_over_e5;
  
  this->current_P_NL = this->P_NL(eta);
  double var_NL_R = variance_of_matter_within_R_NL(R)*var_NL_rescale;
  
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

vector<vector<double> > FlatInhomogeneousUniverseLCDM::return_LOS_integrated_phi_of_lambda(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values){
  
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
  vector<vector<double> > phi_prime_prime_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
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
    dummy_data = compute_phi_tilde_of_lambda_2D(eta, R, f_NL, var_NL_rescale);
    // ISSUE: this order of indeces in confusing.
    y_values[t] = dummy_data[2];
    phi_values[t] = dummy_data[3];
    phi_prime_values[t] = dummy_data[1];
    phi_prime_prime_values[t] = dummy_data[7];
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
    y_values[t] = vector<double>(&y_values[t][0], &y_values[t][lambda_index_y_max]+1);
    phi_values[t] = vector<double>(&phi_values[t][0], &phi_values[t][lambda_index_y_max]+1);
    phi_prime_values[t] = vector<double>(&phi_prime_values[t][0], &phi_prime_values[t][lambda_index_y_max]+1);
    phi_prime_prime_values[t] = vector<double>(&phi_prime_prime_values[t][0], &phi_prime_prime_values[t][lambda_index_y_max]+1);
    
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
      projected_phi_data[3][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &phi_prime_prime_values[t], constants::order_of_interpolation)*pow(kernel_values[i], 2)*dw;
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

double FlatInhomogeneousUniverseLCDM::return_LOS_integrated_variance(double theta, vector<double> w_values, vector<double> kernel_values, double var_NL_rescale){
  
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
  
  return var_projected*var_NL_rescale;
  
}


/*******************************************************************************************************************************************************
 * return_LOS_integrated_skewness
 * Description:
 * - set phi(delta, eta) and lambda(delta, eta) on a 2D grid. This grid is used when computing the LOS-projected CGF in Limber approximation
 * Arguments:
 * - double theta: angular top-hat radius with which LOS-integrated field is smoothed
 * 
*******************************************************************************************************************************************************/

double FlatInhomogeneousUniverseLCDM::return_LOS_integrated_skewness(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values){
  
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
    vNL = variance_of_matter_within_R_NL_2D(R)*var_NL_rescale;
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

double FlatInhomogeneousUniverseLCDM::return_3D_skewness(double z, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale){
  
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
  
  double vNL = variance_of_matter_within_R_NL(R)*var_NL_rescale;
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
    
  vector<vector<double> > phi_data = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, var_NL_rescale, w_values, kernel_values);
  
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
  vector<double> phi_prime_prime_values(n_lambda, 0.0);
  vector<double> tau_values(n_lambda, 0.0);
  
  for(int i = 0; i < tau_values.size(); i++){
    lambda_values[i] = phi_data[0][i];
    phi_values[i] = phi_data[1][i];
    phi_prime_prime_values[i] = phi_data[3][i];
    delta_NL_values[i] = phi_data[2][i];
    if(phi_data[2][i]*phi_data[0][i] > phi_data[1][i]){
      tau_values[i] = sqrt(2.0*(phi_data[2][i]*phi_data[0][i] - phi_data[1][i]));
      if(lambda_values[i] < 0.0) tau_values[i] *= -1.0;
    }
  }
  
  
  
  int n_delta = constants::N_delta_values_for_PDFs;
  double var = interpolate_neville_aitken(0.0, &lambda_values, &phi_prime_prime_values, constants::order_of_interpolation);
  double std = sqrt(var);
  double skew = interpolate_neville_aitken_derivative(0.0, &lambda_values, &phi_prime_prime_values, constants::order_of_interpolation);
  double lognormal_shift = lognormal_tools::get_delta0(var, skew);
  double var_Gauss = log(1.0+var/lognormal_shift/lognormal_shift);
  double mean_Gauss = -0.5*var_Gauss;
  double std_Gauss = sqrt(var_Gauss);
  double delta_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double delta_max = lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0);
  // --> ISSUE: choosing delta_max to be 5*sigma may not be 100% reasonable for very skewed PDFs
  //            Maybe choose it by some quantile in a log-normal PDF that approximates the real PDF?
  
  
  
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
    cout << i << "   ";
    cout << lambda_for_fit[i] << "   ";
    cout << phi_for_fit[i] << "   ";
    cout << phi_prime_for_fit[i] << "\n";
  }
  cout << var << "   ";
  cout << skew << '\n';
  cout << delta_min << '\n';
  cout << delta_max << '\n';
  
  
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
  //for(int i = 0; i < coefficients_phi_of_tau_prime.size()-1; i++) coefficients_phi_of_tau[i+1] = coefficients_phi_of_tau_prime[i]/double(i+1);
  //for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++) coefficients_phi_of_tau_prime[i] = coefficients_phi_of_tau[i+1]*double(i+1);
  
  /*
  for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++){
    cout << i << "   ";
    cout << coefficients_phi_of_tau[i] << "   ";
    cout << coefficients_lambda_of_tau_prime[i] << "\n";
  }
  
  vector<double> coeffs_phi_of_lambda(coefficients_lambda_of_tau.size(), 0.0);
  vector<double> coeffs_phi_prime_of_lambda(coefficients_lambda_of_tau.size(), 0.0);
  //this->return_LOS_integrated_polynomial_coefficients(theta, f_NL, var_NL_rescale, w_values, kernel_values, &coeffs_phi_of_lambda, &coeffs_phi_prime_of_lambda);
  
  vector<vector<double> > Bell_matrix(0, vector<double>(0, 0.0));
  vector<vector<double> > inverse_Bell_matrix(0, vector<double>(0, 0.0));
  vector<double> fact = factoria(coefficients_lambda_of_tau.size());
  return_Bell_matrix(&Bell_matrix, &inverse_Bell_matrix, &coefficients_lambda_of_tau);
  
  for(int i = 0; i < coefficients_lambda_of_tau.size(); i++){
    coeffs_phi_of_lambda[i] = 0.0;
    for(int j = 0; j <= i; j++){
      coeffs_phi_of_lambda[i] += inverse_Bell_matrix[i][j]*coefficients_phi_of_tau[j]*fact[j];
    }
    coeffs_phi_of_lambda[i] /= fact[i];
  }
  
  coeffs_phi_of_lambda[0] = 0.0;
  coeffs_phi_of_lambda[1] = 0.0;
  coeffs_phi_of_lambda[2] = 0.5*var;
  coeffs_phi_of_lambda[3] = skew/6.0;
  
  for(int i = 0; i < coefficients_lambda_of_tau.size(); i++){
    coefficients_phi_of_tau[i] = 0.0;
    coefficients_phi_of_tau_prime[i] = 0.0;
    for(int j = 0; j < coefficients_lambda_of_tau.size(); j++){
      coefficients_phi_of_tau[i] += Bell_matrix[i][j]*coeffs_phi_of_lambda[j]*fact[j];
    }
    coefficients_phi_of_tau[i] /= fact[i];
  }
  
  
  //for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++) coefficients_phi_of_tau[i+1] = coefficients_phi_of_tau[i]/double(i+1);
  for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++) coefficients_phi_of_tau_prime[i] = coefficients_phi_of_tau[i+1]*double(i+1);
  
  
  for(int i = 0; i < coefficients_lambda_of_tau.size()-1; i++){
    cout << i << "   ";
    cout << coeffs_phi_of_lambda[i] << "   ";
    cout << coefficients_phi_of_tau[i] << "   ";
    cout << coefficients_lambda_of_tau_prime[i] << "\n";
  }
  
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
    cout << lambda << "   ";
    cout << tau << "   ";
    cout << exponent << "   ";
    cout << dr << "   ";
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
 * FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF_moment_version_and_saddle_point
 * 
 * Returns PDF of line-of-sight projected matter density contrast with projection kernel
 * specified by the input arrays w_values and kernel_values. w_values should store the
 * lower-redshift edge of each histogram bin and n_of_w_values should contain the redshift
 * histrogram (the histogram doesn't need to be normalised, since normalisation is enforced
 * later on in the code). First column of the returned array is \delta smoothed with a
 * spherical top-hat of R_in_Mpc_over_h, while 2nd column is p(\delta).
 * 
 * However: this function solves the full inverse Laplace transform (for going from CGF to
 * PDF) with a saddle point approximation. And in computed uses the Taylor coefficient of the
 * LOS-integrated CGF instead of the full CGF.
 * 
 */



vector<vector<double> > FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF_saddle_point(vector<double> w_values, vector<double> kernel_values, double theta, double f_NL, double var_NL_rescale){
  
  cout << "Computing projected phi_data:\n";
  cout.flush();
  
  
  vector<vector<double> > phi_data = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, var_NL_rescale, w_values, kernel_values);
  
  /*
   * Determine critical point, where phi(lambda) splits into two branches on the complex plane. 
   * 
   */
  int n_lambda = 0;
  double lambda_c = phi_data[0][0];
  double delta_c = phi_data[2][0];
  for(int i = 1; i < phi_data[0].size(); i++){
    if(phi_data[0][i-1] < phi_data[0][i]){
      n_lambda = i+1;
      lambda_c = phi_data[0][i];
      delta_c = phi_data[2][i];
    }
    else{
      i = 2*phi_data[0].size();
    }
  }
  
  /*
   * Extract phi_data up to the critical point.
   * 
   */
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> phi_values(n_lambda, 0.0);
  vector<double> phi_prime_values(n_lambda, 0.0);
  vector<double> phi_prime_prime_values(n_lambda, 0.0);
  
  for(int i = 0; i < n_lambda; i++){
    lambda_values[i] = phi_data[0][i];
    phi_values[i] = phi_data[1][i];
    phi_prime_values[i] = phi_data[2][i];
    phi_prime_prime_values[i] = phi_data[3][i];
  }
  
  int n_delta = constants::N_delta_values_for_PDFs;
  double var = interpolate_neville_aitken(0.0, &lambda_values, &phi_prime_prime_values, constants::order_of_interpolation);
  double skew = interpolate_neville_aitken_derivative(0.0, &lambda_values, &phi_prime_prime_values, constants::order_of_interpolation);
  double lognormal_shift = lognormal_tools::get_delta0(var, skew);
  double var_Gauss = log(1.0+var/lognormal_shift/lognormal_shift);
  double mean_Gauss = -0.5*var_Gauss;
  double std_Gauss = sqrt(var_Gauss);
  double delta_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double delta_max = lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0);
  /*
   * ISSUE: still avoiding delta=delta_max here, but just for conformity since this was 
   * necessary in compute_LOS_projected_PDF.
   */
  double ddelta = (delta_max-delta_min)/double(n_delta);
  double delta, lambda, phi, phi_prime_prime, phi_prime_prime_prime;
  double var_lambda;
  
  cout << "delta_min = " << delta_min << '\n';
  cout << "delta_max = " << delta_max << '\n';
  
  vector<vector<double> > PDF_data(2, vector<double>(n_delta));
    
  
  /* 
   * Saddle point approximation:
   * 
   * p(\delta) = \int \frac{\mathrm{d}\lambda}{2\pi} \exp(-[i\lambda+\lambda_c]\delta + \phi(i\lambda+\lambda_c))
   * \approx \exp(-\lambda_c\delta + \phi(\lambda_c) x
   * x \int \frac{\mathrm{d}\lambda}{2\pi}  \exp(-i\lambda\delta + i\lambda\phi'(\lambda_c) - 0.5\lambda^2\phi''(\lambda_c))
   * = \exp(-\lambda_c\delta + \phi(\lambda_c) x
   * x \int \frac{\mathrm{d}\lambda}{2\pi}  \exp(-\lambda^2\phi''(\lambda_c))  ; \phi''(\lambda_c) == 1/var_\lambda
   * = \exp(-\lambda_c\delta + \phi(\lambda_c) \sqrt{\frac{1}{2\pi}}
   * x \int \frac{\mathrm{d}\lambda}{\sqrt{2\pi}}  \exp(-0.5\lambda^2\phi''(\lambda_c))
   * = \exp(-\lambda_c\delta + \phi(\lambda_c) / \sqrt{2\pi \phi''(\lambda_c)}
   * 
   */
  cout << "Computing 2D PDF:\n";
  for(int i = 0; i < n_delta; i++){
    delta = delta_min + double(i)*ddelta;
    lambda = interpolate_neville_aitken(delta, &phi_prime_values, &lambda_values, constants::order_of_interpolation);
    phi = interpolate_neville_aitken(delta, &phi_prime_values, &phi_values, constants::order_of_interpolation);
    phi_prime_prime = interpolate_neville_aitken(delta, &phi_prime_values, &phi_prime_prime_values, constants::order_of_interpolation);
    phi_prime_prime_prime = interpolate_neville_aitken_derivative(lambda, &lambda_values, &phi_prime_prime_values, constants::order_of_interpolation);
    PDF_data[0][i] = delta;
    PDF_data[1][i] = exp(-lambda*delta + phi)/sqrt(constants::pi2*phi_prime_prime);
    var_lambda = 1.0/phi_prime_prime;
    //PDF_data[1][i] *= 1.0 + 0.5*pow(phi_prime_prime_prime/6.0, 2)*15.0*pow(var_lambda, 3);
    cout << PDF_data[0][i] << "   ";
    cout << lambda << "   ";
    cout << phi << "   ";
    cout << phi_prime_prime << "   ";
    cout << PDF_data[1][i] << "   ";
    cout << PDF_data[1][i]*(1.0 - 0.5*pow(phi_prime_prime_prime/6.0, 2)*15.0*pow(var_lambda, 3) + 1.0/12.0*pow(phi_prime_prime_prime/6.0, 4)*11.*135135.*pow(var_lambda, 6)) << "\n";
  }
  cout << "Done.\n";
  
  
  
  double var_LOS = this->return_LOS_integrated_variance(theta, w_values, kernel_values, var_NL_rescale);
  cout << "var_LOS = " << var_LOS <<  '\n';
  
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

void FlatInhomogeneousUniverseLCDM::compute_polynomial_coefficients_from_CGF(double e, double R, double f_NL, double var_NL_rescale, vector<double> lambda_values, vector<double> tau_values, vector<double> phi_values, vector<double> phi_prime_values, vector<double> *coeffs_phi_of_lambda, vector<double> *coeffs_phi_prime_of_lambda){
  
  
  
  /*
   * Beginning of old code
   * 
   * 
   */
  
  vector<double> delta_values_temp(0, 0.0);
  vector<double> F_values(0, 0.0);
  vector<double> F_prime_values(0, 0.0);

  this->return_delta_NL_of_delta_L_and_dF_ddelta_2D(e, &delta_values_temp, &F_values, &F_prime_values);
  int n = delta_values_temp.size();  
  double D_growth = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);

  vector<double> tau_prime(n, 0.0);
  vector<double> tau_values_temp(n, 0.0);
  vector<double> y_values(n, 0.0);
  vector<double> phi_values_temp(n, 0.0);
  
  double sigma_R1_squared;
  double dsigma_R1_squared_dlogR;
  double log_R1 = log(R);
  
  for(int i = 0; i < n; i++){
    sigma_R1_squared = interpolate_neville_aitken(log_R1 + 0.5*log(1.0+F_values[i]) - log(c_over_e5), &this->log_top_hat_cylinder_radii, &this->top_hat_cylinder_variances, constants::order_of_interpolation);  
    dsigma_R1_squared_dlogR = interpolate_neville_aitken_derivative(log_R1 + 0.5*log(1.0+F_values[i]) - log(c_over_e5), &this->log_top_hat_cylinder_radii, &this->top_hat_cylinder_variances, constants::order_of_interpolation);  
    delta_values_temp[i] *= D_growth;
    F_prime_values[i] /= D_growth;
    sigma_R1_squared *= D_growth*D_growth;
    dsigma_R1_squared_dlogR *= D_growth*D_growth;
    
    
    tau_values_temp[i] = delta_values_temp[i]/sqrt(sigma_R1_squared);
    
    // Tripple checked:
    tau_prime[i] = (1.0 - 0.25*dsigma_R1_squared_dlogR*F_prime_values[i]/(1.0+F_values[i])*delta_values_temp[i]/sigma_R1_squared)/sqrt(sigma_R1_squared);
    y_values[i] = tau_values_temp[i]*tau_prime[i]/F_prime_values[i];
    phi_values_temp[i] = y_values[i]*F_values[i] - 0.5*pow(tau_values_temp[i], 2.0);
    
  }
  
  
  vector<vector<double> > CGF_data = this->compute_phi_tilde_of_lambda_2D(e, R, f_NL, var_NL_rescale);
  tau_values_temp = CGF_data[0];
  F_values = CGF_data[1];
  y_values = CGF_data[2];
  phi_values_temp = CGF_data[3];
  for(int d = 0; d < tau_values_temp.size(); d++){
    tau_values_temp[d] /= sqrt(CGF_data[4][d]);
  }
  
  
  
  
  int i_min = find_index(0.0, &tau_values_temp)-3;
  int i_max = find_index(0.0, &tau_values_temp)+3;
  while(y_values[i_min-1] < y_values[i_min] && i_min > 0) i_min--;
  while(y_values[i_max+1] > y_values[i_max]  && i_max < n-1) i_max++;

  int n_reduced = i_max - i_min + 1;
  vector<double> tau_values_reduced(n_reduced, 0.0);
  vector<double> y_values_reduced(n_reduced, 0.0);
  vector<double> delta_values_reduced(n_reduced, 0.0);
  vector<double> F_values_reduced(n_reduced, 0.0);

  // This gives the rescaled gen. function phi(y) = Sum <d^n>_{c, NL} y^n.     
  for(int i = 0; i < n_reduced; i++){
    tau_values_reduced[i] = tau_values_temp[i+i_min];
    y_values_reduced[i] = y_values[i+i_min];
    delta_values_reduced[i] = delta_values_temp[i+i_min];
    F_values_reduced[i] = F_values[i+i_min];
  }
  
  
  int n_aux = 30;
  n_aux = 3*constants::generating_function_coeff_order + 1;
  vector<double> tau_aux(n_aux, 0.0);
  vector<double> y_aux(n_aux, 0.0);
  vector<double> G(n_aux, 0.0);
      
  double tau_max = tau_values_temp[i_max];
  double tau_min = tau_values_temp[i_min];
  double dtau = -tau_min/double(n_aux/2);
  for(int i = 0; i < n_aux/2; i++){
    tau_aux[i] = tau_min + double(i)*dtau;
  }
  tau_aux[n_aux/2] = 0.0;
  dtau = tau_max/double(n_aux - n_aux/2 - 1);
  for(int i = n_aux/2 + 1 ; i < n_aux; i++){
    tau_aux[i] = double(i - n_aux/2)*dtau;
  }
  for(int i = 0; i < n_aux; i++){
    y_aux[i] = interpolate_Newton(tau_aux[i] , &tau_values_reduced, &y_values_reduced, constants::order_of_interpolation);
  }

  vector<double> aux_coeffs_y = return_coefficients(&tau_aux, &y_aux, constants::generating_function_coeff_order);
  int polynomial_order = aux_coeffs_y.size()-1;
  
  vector<double> aux_coeffs_G(polynomial_order+1, 0.0);
  vector<double> fact = factoria(polynomial_order+1);
  vector<vector<double> > Bell_matrix(0, vector<double>(0, 0.0));
  vector<vector<double> > inverse_Bell_matrix(0, vector<double>(0, 0.0));
  return_Bell_matrix(&Bell_matrix, &inverse_Bell_matrix, &aux_coeffs_y);
  
  
  
  (*coeffs_phi_prime_of_lambda) = vector<double>(polynomial_order+1, 0.0);
  (*coeffs_phi_of_lambda) = vector<double>(polynomial_order+1, 0.0);

  
  for(int i = 0; i < n_aux; i++){
    G[i] = interpolate_Newton(tau_aux[i] , &tau_values_reduced, &F_values_reduced, constants::order_of_interpolation);
  }
  
  aux_coeffs_G = return_coefficients(&tau_aux, &G, polynomial_order);
  
  for(int i = 0; i < polynomial_order+1; i++){
    for(int j = 0; j <= i; j++){
      (*coeffs_phi_prime_of_lambda)[i] += inverse_Bell_matrix[i][j]*aux_coeffs_G[j]*fact[j];
    }
    (*coeffs_phi_prime_of_lambda)[i] /= fact[i];
  }
  
  for(int i = polynomial_order; i > 0; i--){
    (*coeffs_phi_of_lambda)[i] = (*coeffs_phi_prime_of_lambda)[i-1]/double(i);
  }
  
  this->current_P_L = this->P_L(e);
  double vL = this->variance_of_matter_within_R_NL_2D(R);
  double vL_2 = dsigma_R1_squared_dlogR = interpolate_neville_aitken_derivative(0.0, &y_values_reduced, &F_values_reduced, constants::order_of_interpolation);
  for(int i = 0; i < polynomial_order+1; i++){
    cout << vL << " , ";
    cout << vL_2 << " , ";
    cout << (*coeffs_phi_prime_of_lambda)[i] << " , ";
    cout << (*coeffs_phi_of_lambda)[i] << "\n";
  }
  cout << '\n';
  
  
  /*
   * Ending of old code
   * 
   * 
   */
  
  /*
  
  int n_lambda = lambda_values.size();
  
  int tau_c = tau_values[0];
  for(int i = 0; i < n_lambda-1; i++){
    if(lambda_values[i] < lambda_values[i+1]){
      tau_c = tau_values[i];
    }
  }
  
  n_lambda = y_values.size();
  
  int i_min = find_index(0.0, &tau_values_temp)-3;
  int i_max = find_index(0.0, &tau_values_temp)+3;
  while(y_values[i_min-1] < y_values[i_min] && i_min > 0) i_min--;
  while(y_values[i_max+1] > y_values[i_max]  && i_max < n_lambda-1) i_max++;
  
  int n_reduced = i_max - i_min + 1;
  vector<double> tau_values_reduced(n_reduced, 0.0);
  vector<double> y_values_reduced(n_reduced, 0.0);
  vector<double> F_values_reduced(n_reduced, 0.0);
  
  for(int i = 0; i < n_reduced; i++){
    tau_values_reduced[i] = tau_values_temp[i+i_min];
    y_values_reduced[i] = y_values[i+i_min];
    F_values_reduced[i] = F_values[i+i_min];
  }
  
  
  //double tau_max = 0.9*tau_c;
  //double tau_min = 0.9*tau_values[0];
  double tau_max = tau_values_temp[i_max];
  double tau_min = tau_values_temp[i_min];
  int n_tau = 3*constants::generating_function_coeff_order + 1; // has to be odd number in order to include tau=0 exactly.
  // ISSUE: this shouldn't be hardcoded!
  
  vector<double> tau_for_fit(n_tau,0.0);
  vector<double> lambda_for_fit(n_tau,0.0);
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
    lambda_for_fit[i] = interpolate_Newton(tau_for_fit[i], &tau_values_reduced, &y_values_reduced, constants::order_of_interpolation);
    phi_prime_for_fit[i] = interpolate_Newton(tau_for_fit[i], &tau_values_reduced, &F_values_reduced, constants::order_of_interpolation);
    cout << i << "   ";
    cout << tau_for_fit[i] << "   ";
    cout << lambda_for_fit[i] << "   ";
    cout << phi_prime_for_fit[i] << "\n";
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
    //cout << i << "! = ";
    //cout << fact[i] << "\n";
    for(int j = 0; j <= i; j++){
      (*coeffs_phi_prime_of_lambda)[i] += inverse_Bell_matrix[i][j]*coeffs_phi_prime_of_tau[j]*fact[j];
    }
    (*coeffs_phi_prime_of_lambda)[i] /= fact[i];
  }
  
  
  for(int i = N_coeffs-1; i > 0; i--){
    (*coeffs_phi_of_lambda)[i] = (*coeffs_phi_prime_of_lambda)[i-1]/double(i);
  }
  
  
  for(int i = 0; i < N_coeffs; i++){
    cout << coeffs_lambda_of_tau[i] << " , ";
    cout << coeffs_phi_prime_of_tau[i] << " , ";
    cout << (*coeffs_phi_prime_of_lambda)[i] << " , ";
    cout << (*coeffs_phi_of_lambda)[i] << "\n";
  }
  
  */
  
}





/*******************************************************************************************************************************************************
 * return_LOS_integrated_polynomial_coefficients
 * 
 * Returns polynomial coefficients of the cumulant generating function, projected along the line-of-sight kernel kernel_values.
 * NOTE: this function does not check whether the binning of w_values is too coarse for accurate integration.
 *       (The binning is instead set and controlled by ProjectedGalaxySample::set_n_of_w_data)
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::return_LOS_integrated_polynomial_coefficients(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values, vector<double> *coeffs_phi_of_lambda, vector<double> *coeffs_phi_prime_of_lambda){
  
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
  vector<double> lambda_values;
  vector<double> tau_values;
  vector<double> phi_values;
  vector<double> phi_prime_values;
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time_of_nonzero_kernel; t++){
    int i = indeces_of_nonzero_kernel[t];
    w = w_values_bin_center[i];
    eta = eta_0-w;
    R = w*theta;
    
    CGF_data = this->compute_phi_tilde_of_lambda_2D(eta, R, f_NL, var_NL_rescale);
    tau_values = CGF_data[0];
    lambda_values = CGF_data[2];
    phi_values = CGF_data[3];
    phi_prime_values = CGF_data[1];
    for(int d = 0; d < tau_values.size(); d++){
      tau_values[d] /= sqrt(CGF_data[4][d]);
    }
    this->compute_polynomial_coefficients_from_CGF(eta, R, f_NL, var_NL_rescale, lambda_values, tau_values, phi_values, phi_prime_values, &coefficients_phi_of_lambda[t], &coefficients_phi_prime_of_lambda[t]);
    coefficients_phi_prime_of_lambda[t][0] = 0.0;
    // ISSUE: you can do this because current non-linear P(k) has been computed in "compute_phi_tilde_of_lambda_2D".
    // But this kind of flow seems very prone to error!
    coefficients_phi_prime_of_lambda[t][1] = variance_of_matter_within_R_NL_2D(R)*var_NL_rescale;
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

void FlatInhomogeneousUniverseLCDM::return_LOS_integrated_polynomial_coefficients_incl_CMB_kappa(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values, vector<vector<double> > *coeffs_phi, vector<vector<double> > *coeffs_dphi_dlambda_delta, vector<vector<double> > *coeffs_dphi_dlambda_kappa){
  
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
    
    CGF_data = this->compute_phi_tilde_of_lambda_2D(eta, R, f_NL, var_NL_rescale);
    tau_values = CGF_data[0];
    for(int d = 0; d < tau_values.size(); d++){
      tau_values[d] /= sqrt(CGF_data[4][d]);
    }
    this->compute_polynomial_coefficients_from_CGF(eta, R, f_NL, var_NL_rescale, CGF_data[2], tau_values, CGF_data[3], CGF_data[1], &coefficients_phi_of_lambda[t], &coefficients_phi_prime_of_lambda[t]);
    coefficients_phi_prime_of_lambda[t][0] = 0.0;
    // ISSUE: you can do this because current non-linear P(k) has been computed in "compute_phi_tilde_of_lambda_2D".
    // But this kind of flow seems very prone to error!
    coefficients_phi_prime_of_lambda[t][1] = variance_of_matter_within_R_NL_2D(R)*var_NL_rescale;
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

/*******************************************************************************************************************************************************
 * return_LOS_integrated_phi_of_lambda_incl_CMB_kappa
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::return_LOS_integrated_phi_of_lambda_incl_CMB_kappa(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values, vector<vector<double> > *phi_data_delta, vector<vector<double> > *phi_data_kappa, vector<vector<double> > *phi_grid, vector<vector<double> > *dphi_dldelta_grid, vector<vector<double> > *dphi_dlkappa_grid, vector<vector<double> > *d2phi_dldelta2_grid, vector<vector<double> > *d2phi_dldelta_dlkappa_grid, vector<vector<double> > *d2phi_dlkappa2_grid, vector<vector<int> > *grid_mask){
  
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
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  vector<vector<double> > y_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_prime_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_prime_prime_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > dummy_data;
  
  double y_min = 0.0;
  double y_max = 0.0;
  vector<double>::iterator y_min_iterator;
  vector<double>::iterator y_max_iterator;
  int time_index_y_max, lambda_index_y_max, lambda_index_y_min;
  
  cout << "HAHAHAHAH 1\n";
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time; t++){
    cout << w << "  ";
    cout.flush();
    w = w_values_bin_center[t];
    eta = eta_0-w;
    R = w*theta;
    cout << w << "  ";
    cout.flush();
    dummy_data = compute_phi_tilde_of_lambda_2D(eta, R, f_NL, var_NL_rescale);
    cout << t << "  ";
    cout.flush();
    // ISSUE: this order of indeces in confusing.
    y_values[t] = dummy_data[2];
    phi_values[t] = dummy_data[3];
    phi_prime_values[t] = dummy_data[1];
    phi_prime_prime_values[t] = dummy_data[7];
    
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
    y_values[t] = vector<double>(&y_values[t][0], &y_values[t][lambda_index_y_max]+1);
    phi_values[t] = vector<double>(&phi_values[t][0], &phi_values[t][lambda_index_y_max]+1);
    phi_prime_values[t] = vector<double>(&phi_prime_values[t][0], &phi_prime_values[t][lambda_index_y_max]+1);
    phi_prime_prime_values[t] = vector<double>(&phi_prime_prime_values[t][0], &phi_prime_prime_values[t][lambda_index_y_max]+1);
    
    cout << n_time << '\n';
    cout.flush();
  }
  
  cout << "computing ranges of lambda_delta and lambda_kappa & resizing phi_grids accordingly\n";
  
  vector<int> indeces_of_nonzero_kernel(0,0);
  vector<double> CMB_lensing_kernel_nonoverlap = CMB_lensing_kernel;
  for(int t = 0; t < n_time; t++){
    if(kernel_values_extended_to_last_scattering[t] > 0.0){
      indeces_of_nonzero_kernel.push_back(t);
      CMB_lensing_kernel_nonoverlap[t] = 0.0;
    }
  }
  
  (*phi_data_delta) = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, var_NL_rescale, w_values, kernel_values);
  (*phi_data_kappa) = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, var_NL_rescale, w_values_extended_to_last_scattering, CMB_lensing_kernel);
  vector<vector<double> > phi_data_kappa_nonoverlap = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, var_NL_rescale, w_values_extended_to_last_scattering, CMB_lensing_kernel_nonoverlap);
  
  
  double var_delta = interpolate_neville_aitken(0.0, &(*phi_data_delta)[0], &(*phi_data_delta)[3], constants::order_of_interpolation);
  double std_delta = sqrt(var_delta);
  double skew = interpolate_neville_aitken_derivative(0.0, &(*phi_data_delta)[0], &(*phi_data_delta)[3], constants::order_of_interpolation);
  double lognormal_shift = lognormal_tools::get_delta0(var_delta, skew);
  double var_Gauss = log(1.0+var_delta/lognormal_shift/lognormal_shift);
  double mean_Gauss = -0.5*var_Gauss;
  double std_Gauss = sqrt(var_Gauss);
  double delta_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double delta_max = lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0);
  
  double var_kappa = interpolate_neville_aitken(0.0, &(*phi_data_kappa)[0], &(*phi_data_kappa)[3], constants::order_of_interpolation);
  double std_kappa = sqrt(var_kappa);
  skew = interpolate_neville_aitken_derivative(0.0, &(*phi_data_kappa)[0], &(*phi_data_kappa)[3], constants::order_of_interpolation);
  lognormal_shift = lognormal_tools::get_delta0(var_kappa, skew);
  var_Gauss = log(1.0+var_kappa/lognormal_shift/lognormal_shift);
  mean_Gauss = -0.5*var_Gauss;
  std_Gauss = sqrt(var_Gauss);
  double kappa_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double kappa_max = lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0);
  
  
  int N_ldelta = (*phi_data_delta)[0].size();
  int N_lkappa = (*phi_data_kappa)[0].size();
  (*phi_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*dphi_dldelta_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*dphi_dlkappa_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*d2phi_dldelta2_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*d2phi_dldelta_dlkappa_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*d2phi_dlkappa2_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*grid_mask) = vector<vector<int> >(N_ldelta, vector<int>(N_lkappa, 0));
  
  int t, N_nonzero_kernel = indeces_of_nonzero_kernel.size();
  double y, ldelta, lkappa;
  double kernel_delta, kernel_kappa;
  double dw_times_kernel_delta;
  double dw_times_kernel_kappa;
  
  double phi;
  double dphi_dldelta;
  double dphi_dlkappa;
  double d2phi_dldelta2;
  double d2phi_dldelta_dlkappa;
  double d2phi_dlkappa2;
  
  double phi_prime_integrand;
  double phi_prime_prime_integrand;
  vector<double> * pointer_to_y_values;
  
  for(int d = 0; d < N_ldelta; d++){
    cout << d << '\n';
    ldelta = (*phi_data_delta)[0][d];
      for(int k = 0; k < N_lkappa; k++){
        lkappa = (*phi_data_kappa)[0][k];
        
        dphi_dldelta = 0.0;
        dphi_dlkappa = interpolate_neville_aitken(lkappa, &phi_data_kappa_nonoverlap[0], &phi_data_kappa_nonoverlap[2], constants::order_of_interpolation);
        
          for(int i = 0; i < N_nonzero_kernel; i++){
            t = indeces_of_nonzero_kernel[i];
            dw = w_values_extended_to_last_scattering[t+1]-w_values_extended_to_last_scattering[t];
            kernel_delta = kernel_values_extended_to_last_scattering[t];
            dw_times_kernel_delta = kernel_delta*dw;
            kernel_kappa = CMB_lensing_kernel[t];
            dw_times_kernel_kappa = kernel_kappa*dw;
            y = kernel_delta*ldelta;
            y += kernel_kappa*lkappa;
            
            phi_prime_integrand = interpolate_neville_aitken(y, &y_values[t], &phi_prime_values[t], constants::order_of_interpolation);
            dphi_dldelta += dw_times_kernel_delta*phi_prime_integrand;
            dphi_dlkappa += dw_times_kernel_kappa*phi_prime_integrand;
            
          }
          
          (*dphi_dldelta_grid)[d][k] = dphi_dldelta;
          (*dphi_dlkappa_grid)[d][k] = dphi_dlkappa;
          
        }
  }
  
  
  for(int d = 0; d < N_ldelta; d++){
    cout << d << '\n';
    ldelta = (*phi_data_delta)[0][d];
      for(int k = 0; k < N_lkappa; k++){
        lkappa = (*phi_data_kappa)[0][k];
        
        phi = interpolate_neville_aitken(lkappa, &phi_data_kappa_nonoverlap[0], &phi_data_kappa_nonoverlap[1], constants::order_of_interpolation);;
        dphi_dldelta = (*dphi_dldelta_grid)[d][k];
        dphi_dlkappa = (*dphi_dlkappa_grid)[d][k];
        d2phi_dldelta2 = 0.0;
        d2phi_dldelta_dlkappa = 0.0;
        d2phi_dlkappa2 = interpolate_neville_aitken(lkappa, &phi_data_kappa_nonoverlap[0], &phi_data_kappa_nonoverlap[3], constants::order_of_interpolation);;
        
        if(dphi_dldelta > delta_min && dphi_dldelta < delta_max && dphi_dlkappa > kappa_min && dphi_dlkappa < kappa_max){
        
          for(int i = 0; i < N_nonzero_kernel; i++){
            t = indeces_of_nonzero_kernel[i];
            dw = w_values_extended_to_last_scattering[t+1]-w_values_extended_to_last_scattering[t];
            kernel_delta = kernel_values_extended_to_last_scattering[t];
            dw_times_kernel_delta = kernel_delta*dw;
            kernel_kappa = CMB_lensing_kernel[t];
            dw_times_kernel_kappa = kernel_kappa*dw;
            y = kernel_delta*ldelta;
            y += kernel_kappa*lkappa;
            
            pointer_to_y_values = &y_values[t];
            lambda_index_y_max = pointer_to_y_values->size() -1;
            if(y < y_values[t][lambda_index_y_max] && y > y_values[t][0]){
              phi += dw*interpolate_neville_aitken(y, pointer_to_y_values, &phi_values[t], constants::order_of_interpolation);
              phi_prime_prime_integrand = interpolate_neville_aitken(y, pointer_to_y_values, &phi_prime_prime_values[t], constants::order_of_interpolation);
              d2phi_dldelta2 += dw_times_kernel_delta*kernel_delta*phi_prime_prime_integrand;
              d2phi_dldelta_dlkappa += dw_times_kernel_delta*kernel_kappa*phi_prime_prime_integrand;
              d2phi_dlkappa2 += dw_times_kernel_kappa*kernel_kappa*phi_prime_prime_integrand;
              (*grid_mask)[d][k] = 1;
            }
          }
          
          (*phi_grid)[d][k] = phi;
          (*d2phi_dldelta2_grid)[d][k] = d2phi_dldelta2;
          (*d2phi_dldelta_dlkappa_grid)[d][k] = d2phi_dldelta_dlkappa;
          (*d2phi_dlkappa2_grid)[d][k] = d2phi_dlkappa2;
          
        }
      }
  }
  
}




/*
 * FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point
 * 
 * 
 * 
 */



void FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point(vector<double> w_values, vector<double> kernel_values, double theta, double f_NL, double var_NL_rescale, vector<vector<double> > *delta_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid){
  
  cout << "Computing projected phi_data:\n";
  cout.flush();
  
  vector<vector<double> > phi_data_delta;
  vector<vector<double> > phi_data_kappa;
  vector<vector<double> > phi_grid;
  vector<vector<double> > dphi_dldelta_grid;
  vector<vector<double> > dphi_dlkappa_grid;
  vector<vector<double> > d2phi_dldelta2_grid;
  vector<vector<double> > d2phi_dldelta_dlkappa_grid;
  vector<vector<double> > d2phi_dlkappa2_grid;
  vector<vector<int> > grid_mask;
  
  this->return_LOS_integrated_phi_of_lambda_incl_CMB_kappa(theta, f_NL, var_NL_rescale, w_values, kernel_values, &phi_data_delta, &phi_data_kappa, &phi_grid, &dphi_dldelta_grid, &dphi_dlkappa_grid, &d2phi_dldelta2_grid, &d2phi_dldelta_dlkappa_grid, &d2phi_dlkappa2_grid, &grid_mask);
  
  
  double var_delta = interpolate_neville_aitken(0.0, &phi_data_delta[0], &phi_data_delta[3], constants::order_of_interpolation);
  double std_delta = sqrt(var_delta);
  double skew = interpolate_neville_aitken_derivative(0.0, &phi_data_delta[0], &phi_data_delta[3], constants::order_of_interpolation);
  double lognormal_shift = lognormal_tools::get_delta0(var_delta, skew);
  double var_Gauss = log(1.0+var_delta/lognormal_shift/lognormal_shift);
  double mean_Gauss = -0.5*var_Gauss;
  double std_Gauss = sqrt(var_Gauss);
  double delta_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double delta_max = lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0);
  
  double var_kappa = interpolate_neville_aitken(0.0, &phi_data_kappa[0], &phi_data_kappa[3], constants::order_of_interpolation);
  double std_kappa = sqrt(var_kappa);
  skew = interpolate_neville_aitken_derivative(0.0, &phi_data_kappa[0], &phi_data_kappa[3], constants::order_of_interpolation);
  lognormal_shift = lognormal_tools::get_delta0(var_kappa, skew);
  var_Gauss = log(1.0+var_kappa/lognormal_shift/lognormal_shift);
  mean_Gauss = -0.5*var_Gauss;
  std_Gauss = sqrt(var_Gauss);
  double kappa_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double kappa_max = lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0);
  
  cout << "delta_min = " << delta_min << '\n';
  cout << "delta_max = " << delta_max << '\n';
  cout << "kappa_min = " << kappa_min << '\n';
  cout << "kappa_max = " << kappa_max << '\n';
  
  int n_delta = constants::N_delta_values_for_PDFs;
  int n_kappa = constants::N_delta_values_for_PDFs;
  /*
   * ISSUE: still avoiding delta=delta_max here, but just for conformity since this was 
   * necessary in compute_LOS_projected_PDF.
   */
  double ddelta = (delta_max-delta_min)/double(n_delta);
  double dkappa = (kappa_max-kappa_min)/double(n_kappa);
  
  vector<double> delta_values(n_delta, 0.0);
  vector<double> kappa_values(n_kappa, 0.0);
  for(int i = 0; i < n_delta; i++){
    delta_values[i] = delta_min + double(i)*ddelta;
  }
  for(int i = 0; i < n_kappa; i++){
    kappa_values[i] = kappa_min + double(i)*dkappa;
  }
  
  
  double y, delta, kappa, lambda_delta, lambda_kappa;
  double phi;
  double dphi_dldelta;
  double dphi_dlkappa;
  double d2phi_dldelta2;
  double d2phi_dldelta_dlkappa;
  double d2phi_dlkappa_dldelta;
  double d2phi_dlkappa2;
  double determinant;
  vector<vector<double> > inverse_Jacobian(2, vector<double>(2, 0.0));
  
  
  (*delta_grid) = vector<vector<double> >(n_delta, vector<double>(n_kappa, 0.0));
  (*kappa_grid) = vector<vector<double> >(n_delta, vector<double>(n_kappa, 0.0));
  (*PDF_grid) = vector<vector<double> >(n_delta, vector<double>(n_kappa, 0.0));
  
  int steps;
  
  for(int d = 0; d < n_delta; d++){
    delta = delta_values[d];
    cout << delta << "   ";
    for(int k = 0; k < n_kappa; k++){
      kappa = kappa_values[k];
      (*delta_grid)[d][k] = delta;
      (*kappa_grid)[d][k] = kappa;
      
      lambda_kappa = 0.0;//interpolate_neville_aitken(kappa, &phi_data_kappa[2], &phi_data_kappa[0], constants::order_of_interpolation);
      lambda_delta = 0.0;//interpolate_neville_aitken(delta, &phi_data_delta[2], &phi_data_delta[0], constants::order_of_interpolation);
      
      dphi_dldelta = interpolate_neville_aitken_grid(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &dphi_dldelta_grid, constants::order_of_interpolation, constants::order_of_interpolation);
      dphi_dlkappa = interpolate_neville_aitken_grid(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &dphi_dlkappa_grid, constants::order_of_interpolation, constants::order_of_interpolation);
      steps = 0;
      // ISSUE: these precision criteria should be set in constants.h in the "constants" namespace.
      while(steps < 1000 && (abs(dphi_dldelta - delta) > 0.1*ddelta || abs(dphi_dlkappa - kappa) > 0.1*dkappa)){
        steps++;
        d2phi_dldelta2 = interpolate_neville_aitken_dgrid_dt(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &dphi_dldelta_grid, constants::order_of_interpolation, constants::order_of_interpolation);
        d2phi_dldelta_dlkappa = interpolate_neville_aitken_dgrid_dx(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &dphi_dldelta_grid, constants::order_of_interpolation, constants::order_of_interpolation);
        d2phi_dlkappa_dldelta = interpolate_neville_aitken_dgrid_dt(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &dphi_dlkappa_grid, constants::order_of_interpolation, constants::order_of_interpolation);
        d2phi_dlkappa2 = interpolate_neville_aitken_dgrid_dx(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &dphi_dlkappa_grid, constants::order_of_interpolation, constants::order_of_interpolation);
        determinant = d2phi_dldelta2*d2phi_dlkappa2;
        determinant -= d2phi_dldelta_dlkappa*d2phi_dlkappa_dldelta;
        
        inverse_Jacobian[0][0] = d2phi_dlkappa2/determinant;
        inverse_Jacobian[0][1] = -d2phi_dldelta_dlkappa/determinant;
        inverse_Jacobian[1][0] = -d2phi_dlkappa_dldelta/determinant;
        inverse_Jacobian[1][1] = d2phi_dldelta2/determinant;
        
        lambda_delta -= inverse_Jacobian[0][0]*(dphi_dldelta-delta) + inverse_Jacobian[0][1]*(dphi_dlkappa-kappa);
        lambda_kappa -= inverse_Jacobian[1][0]*(dphi_dldelta-delta) + inverse_Jacobian[1][1]*(dphi_dlkappa-kappa);
        dphi_dldelta = interpolate_neville_aitken_grid(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &dphi_dldelta_grid, constants::order_of_interpolation, constants::order_of_interpolation);
        dphi_dlkappa = interpolate_neville_aitken_grid(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &dphi_dlkappa_grid, constants::order_of_interpolation, constants::order_of_interpolation);
        
        // f = f0
        // 
        //(f-f0)/(x-x0) = f'
        //==> x0 = x-(f-f0)/f'
        
        // d_lambda = Jac^{-1} * ((delta,kappa) - (dphi_dldelta, dphi_dlkappa))
        
      }
      phi = interpolate_neville_aitken_grid(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &phi_grid, constants::order_of_interpolation, constants::order_of_interpolation);
      d2phi_dldelta2 = interpolate_neville_aitken_grid(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &d2phi_dldelta2_grid, constants::order_of_interpolation, constants::order_of_interpolation);
      d2phi_dldelta_dlkappa = interpolate_neville_aitken_grid(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &d2phi_dldelta_dlkappa_grid, constants::order_of_interpolation, constants::order_of_interpolation);
      d2phi_dlkappa2 = interpolate_neville_aitken_grid(lambda_delta, lambda_kappa, &phi_data_delta[0], &phi_data_kappa[0], &d2phi_dlkappa2_grid, constants::order_of_interpolation, constants::order_of_interpolation);
      determinant = d2phi_dldelta2*d2phi_dlkappa2;
      determinant -= pow(d2phi_dldelta_dlkappa, 2);
      
      int index_delta = find_index(lambda_delta, &phi_data_delta[0]);
      int index_kappa = find_index(lambda_kappa, &phi_data_kappa[0]);
      if(determinant > 0.0 && steps < 1000 && grid_mask[index_delta][index_kappa] == 1){
        (*PDF_grid)[d][k] = exp(-delta*lambda_delta-kappa*lambda_kappa+phi);
        (*PDF_grid)[d][k] /= constants::pi2*sqrt(determinant);
      }
      
    }
    
    cout << steps << '\n';
  }
  
  
  
  /*
  int N_ldelta = phi_data_delta[0].size();
  int N_lkappa = phi_data_kappa[0].size();
  (*PDF_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  
  for(int d = 0; d < N_ldelta; d++){
    lambda_delta = phi_data_delta[0][d];
    for(int k = 0; k < N_lkappa; k++){
      lambda_kappa = phi_data_kappa[0][k];
      delta = (*delta_grid)[d][k];
      kappa = (*kappa_grid)[d][k];
      if(grid_mask[d][k] == 1 && delta > -0.5 && delta < 1.0 && kappa>-0.4 && kappa < 0.4){
        (*PDF_grid)[d][k] = exp(-delta*lambda_delta-kappa*lambda_kappa+phi_grid[d][k]);
        determinant = d2phi_dldelta2_grid[d][k]*d2phi_dlkappa2_grid[d][k];
        determinant -= pow(d2phi_dldelta_dlkappa_grid[d][k], 2);
        if(determinant>0.0)
          (*PDF_grid)[d][k] /= constants::pi2*sqrt(determinant);
        else
          (*PDF_grid)[d][k] = 0.0;
      }
    }
  }
  */
  
}



