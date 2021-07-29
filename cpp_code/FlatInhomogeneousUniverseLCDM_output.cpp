
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


void FlatInhomogeneousUniverseLCDM::return_delta_NL_of_delta_L_and_dF_ddelta_3D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values, vector<double> *delta_NL_prime_prime_values){
    
  int N_delta = this->delta_values_for_spherical_collapse.size();
  
  vector<double> delta_L_dummy = this->delta_values_for_spherical_collapse;
  vector<double> delta_NL_dummy = vector<double>(N_delta);
  vector<double> delta_NL_prime_dummy = vector<double>(N_delta);
  vector<double> delta_NL_prime_prime_dummy = vector<double>(N_delta);
  
  int N_delta_uncollapsed = 0;
  for(int i = 0; i < N_delta; i++){
    delta_NL_dummy[i] = interpolate_neville_aitken(eta, &this->eta, &this->spherical_collapse_evolution_of_delta[i], constants::order_of_interpolation);
    delta_NL_prime_dummy[i] = interpolate_neville_aitken(eta, &this->eta, &this->spherical_collapse_evolution_of_delta_ddelta[i], constants::order_of_interpolation);
    delta_NL_prime_prime_dummy[i] = interpolate_neville_aitken(eta, &this->eta, &this->spherical_collapse_evolution_of_delta_ddelta2[i], constants::order_of_interpolation);
    if(delta_NL_dummy[i] < constants::max_contrast){
      N_delta_uncollapsed++;
    }
    else{
      i = N_delta;
      break;
    }
  }
  
  N_delta_uncollapsed -= constants::order_of_interpolation;
  
  (*delta_L_values) = vector<double>(N_delta_uncollapsed);
  (*delta_NL_values) = vector<double>(N_delta_uncollapsed);
  (*delta_NL_prime_values) = vector<double>(N_delta_uncollapsed);
  (*delta_NL_prime_prime_values) = vector<double>(N_delta_uncollapsed);
  
  for(int i = 0; i < N_delta_uncollapsed; i++){
    (*delta_L_values)[i] = delta_L_dummy[i];
    (*delta_NL_values)[i] = delta_NL_dummy[i];
    (*delta_NL_prime_values)[i] = delta_NL_prime_dummy[i];
    (*delta_NL_prime_prime_values)[i] = delta_NL_prime_prime_dummy[i];
  }
  
}


void FlatInhomogeneousUniverseLCDM::return_delta_NL_of_delta_L_and_dF_ddelta_2D(double eta, vector<double> *delta_L_values, vector<double> *delta_NL_values, vector<double> *delta_NL_prime_values, vector<double> *delta_NL_prime_prime_values){
  
  int N_delta = this->delta_values_for_cylindrical_collapse.size();
  
  vector<double> delta_L_dummy = this->delta_values_for_cylindrical_collapse;
  vector<double> delta_NL_dummy = vector<double>(N_delta);
  vector<double> delta_NL_prime_dummy = vector<double>(N_delta);
  vector<double> delta_NL_prime_prime_dummy = vector<double>(N_delta);
  
  int N_delta_uncollapsed = 0;
  for(int i = 0; i < N_delta; i++){
    delta_NL_dummy[i] = interpolate_neville_aitken(eta, &this->eta, &this->cylindrical_collapse_evolution_of_delta[i], constants::order_of_interpolation);
    delta_NL_prime_dummy[i] = interpolate_neville_aitken(eta, &this->eta, &this->cylindrical_collapse_evolution_of_delta_ddelta[i], constants::order_of_interpolation);
    delta_NL_prime_prime_dummy[i] = interpolate_neville_aitken(eta, &this->eta, &this->cylindrical_collapse_evolution_of_delta_ddelta2[i], constants::order_of_interpolation);
    if(delta_NL_dummy[i] < constants::max_contrast){
      N_delta_uncollapsed++;
    }
    else{
      i = N_delta;
      break;
    }
  }
  
  N_delta_uncollapsed -= constants::order_of_interpolation;
  
  (*delta_L_values) = vector<double>(N_delta_uncollapsed);
  (*delta_NL_values) = vector<double>(N_delta_uncollapsed);
  (*delta_NL_prime_values) = vector<double>(N_delta_uncollapsed);
  (*delta_NL_prime_prime_values) = vector<double>(N_delta_uncollapsed);
  
  for(int i = 0; i < N_delta_uncollapsed; i++){
    (*delta_L_values)[i] = delta_L_dummy[i];
    (*delta_NL_values)[i] = delta_NL_dummy[i];
    (*delta_NL_prime_values)[i] = delta_NL_prime_dummy[i];
    (*delta_NL_prime_prime_values)[i] = delta_NL_prime_prime_dummy[i];
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

vector<vector<double> > FlatInhomogeneousUniverseLCDM::return_power_spectra(double e){
  
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
  double scale, e, t_val, dF_ddelta, d2F_ddelta2;
  fstream growth_stream;
  
  vector<double> delta_L_values;
  vector<double> delta_NL_values;
  vector<double> delta_NL_prime_values;
  vector<double> delta_NL_prime_prime_values;
  
  remove(file_name.c_str());
  FILE * F = fopen(file_name.c_str(), "w");
  fclose(F);
  growth_stream.open(file_name.c_str());
  
  growth_stream << setprecision(10) << scientific;
  growth_stream << "#Cosmological Parameters: Om_m = " << Om_m << ", Om_l = " << Om_l << ", Om_r = " << Om_r << '\n';
  growth_stream << "#a_max = " << this->a_final << '\n';
  growth_stream << "#eta(t)" << setw(20) << "t(eta)" << setw(20) << "w(eta)" << setw(20) << "z(eta)" << setw(20) << "a(eta)" << setw(20) << "D(eta)" << setw(20) << "D_prime(eta)" << setw(20) << " dF_ddelta_cyl(eta)" << setw(20) << " d2F_ddelta2_cyl(eta)\n";
  for(int i = 0; i < this->return_number_of_time_steps(); i++){
    t_val = this->t[i];
    e = this->eta[i];
    scale = this->a[i];
    this->return_delta_NL_of_delta_L_and_dF_ddelta_2D(e, &delta_L_values, &delta_NL_values, &delta_NL_prime_values, &delta_NL_prime_prime_values);
    dF_ddelta = interpolate_neville_aitken(0.0, &delta_L_values, &delta_NL_prime_values, constants::order_of_interpolation);
    d2F_ddelta2 = interpolate_neville_aitken(0.0, &delta_L_values, &delta_NL_prime_prime_values, constants::order_of_interpolation);
    growth_stream << setw(20) << e;
    growth_stream << setw(20) << t_val;
    growth_stream << setw(20) << this->eta_at_a(1.0)-e;
    growth_stream << setw(20) << 1.0/scale - 1.0;
    growth_stream << setw(20) << scale;
    growth_stream << setw(20) << this->Newtonian_growth_factor_of_delta[i];
    growth_stream << setw(20) << this->Newtonian_growth_factor_of_delta_prime[i];
    growth_stream << setw(20) << dF_ddelta;
    growth_stream << setw(20) << d2F_ddelta2 << '\n';
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
  vector<vector<double> > phi_prime_prime_prime_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
  vector<vector<double> > bias_term_1_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
  vector<vector<double> > bias_term_2_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
  vector<vector<double> > variance_term_1_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
  vector<vector<double> > variance_term_2_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
  vector<vector<double> > variance_term_3_values(n_time_of_nonzero_kernel, vector<double>(n_lambda, 0.0));
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
    y_values[t] = dummy_data[0];
    phi_values[t] = dummy_data[1];
    phi_prime_values[t] = dummy_data[2];
    phi_prime_prime_values[t] = dummy_data[3];
    phi_prime_prime_prime_values[t] = dummy_data[4];
    bias_term_1_values[t] = dummy_data[5];
    bias_term_2_values[t] = dummy_data[6];
    variance_term_1_values[t] = dummy_data[7];
    variance_term_2_values[t] = dummy_data[8];
    variance_term_3_values[t] = dummy_data[9];
    n_lambda = y_values[t].size();
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
    phi_prime_prime_prime_values[t] = vector<double>(&phi_prime_prime_prime_values[t][0], &phi_prime_prime_prime_values[t][lambda_index_y_max]+1);
    bias_term_1_values[t] = vector<double>(&bias_term_1_values[t][0], &bias_term_1_values[t][lambda_index_y_max]+1);
    bias_term_2_values[t] = vector<double>(&bias_term_2_values[t][0], &bias_term_2_values[t][lambda_index_y_max]+1);
    variance_term_1_values[t] = vector<double>(&variance_term_1_values[t][0], &variance_term_1_values[t][lambda_index_y_max]+1);
    variance_term_2_values[t] = vector<double>(&variance_term_2_values[t][0], &variance_term_2_values[t][lambda_index_y_max]+1);
    variance_term_3_values[t] = vector<double>(&variance_term_3_values[t][0], &variance_term_3_values[t][lambda_index_y_max]+1);
    
  }
  
  lambda_index_y_min = 0;
  lambda_index_y_max = y_values[time_index_y_max].size()-1;
  while(y_values[time_index_y_max][lambda_index_y_min] < y_min){
    lambda_index_y_min++;
  }
  
  int n_lambda_1st_branch = 1 + lambda_index_y_max - lambda_index_y_min;
  /*
   * 0 = lambda
   * 1 = phi
   * 2 = phi'
   * 3 = phi''
   * 4 = phi'''
   * 5 = bias_term_1
   * 6 = bias_term_2
   * 7 = variance_term_1
   * 8 = variance_term_2
   * 9 = variance_term_3
   * 
   */
  vector<vector<double> > projected_phi_data(10, vector<double>(n_lambda_1st_branch, 0.0));
  projected_phi_data[0] = vector<double>(&y_values[time_index_y_max][lambda_index_y_min], &y_values[time_index_y_max][lambda_index_y_max]+1);
  
  cout << "projecting CGF\n";
  
  double delta_L, delta_NL;
  double kernel_times_dw;
  double kernelSq_times_dw;
  for(int y = 0; y < n_lambda_1st_branch; y++){
    for(int t = 0; t < n_time_of_nonzero_kernel; t++){
      int i = indeces_of_nonzero_kernel[t];
      dw = w_values[i+1]-w_values[i];
      kernel_times_dw = dw*kernel_values[i];
      kernelSq_times_dw = kernel_times_dw*kernel_values[i];
      projected_phi_data[1][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &phi_values[t], constants::order_of_interpolation)*dw;
      projected_phi_data[2][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &phi_prime_values[t], constants::order_of_interpolation)*kernel_times_dw;
      projected_phi_data[3][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &phi_prime_prime_values[t], constants::order_of_interpolation)*kernelSq_times_dw;
      projected_phi_data[4][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &phi_prime_prime_prime_values[t], constants::order_of_interpolation)*kernel_values[i]*kernelSq_times_dw;
      // Now bias terms:
      // <delta_g | delta_m > = invLaplace[delta_m + b_1_Lagrange*bias_term_1 + b_2_Lagrange*bias_term_2]
      projected_phi_data[5][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &bias_term_1_values[t], constants::order_of_interpolation)*kernel_times_dw;
      projected_phi_data[6][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &bias_term_2_values[t], constants::order_of_interpolation)*kernel_times_dw;
      // Now bias terms:
      // <delta_g^2 | delta_m > = invLaplace[var_term_1 + b_1_Lagrange*var_term_2 + b_1_Lagrange^2*var_term_3]
      // ISSUE: what happens to 'L' in Limber approx?
      projected_phi_data[7][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &variance_term_1_values[t], constants::order_of_interpolation)*kernelSq_times_dw;
      projected_phi_data[8][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &variance_term_2_values[t], constants::order_of_interpolation)*kernelSq_times_dw;
      projected_phi_data[9][y] += interpolate_neville_aitken(projected_phi_data[0][y], &y_values[t], &variance_term_3_values[t], constants::order_of_interpolation)*kernelSq_times_dw;
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
    
  vector<vector<double> > CGF_data = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, var_NL_rescale, w_values, kernel_values);
  
  cout << "Computing CGF_data done.\n";
  
  int n_lambda = 0;
  for(int i = 1; i < CGF_data[0].size(); i++){
    // ISSUE: this may fail for strong PNG, when there is two branch cuts in the CGF
    if(CGF_data[0][i-1] < CGF_data[0][i]) n_lambda = i+1;
  }
  
  vector<double> lambda_values(n_lambda, 0.0);
  vector<double> d2phi_dlambda2(n_lambda, 0.0);
  for(int i = 0; i < n_lambda; i++){
    lambda_values[i] = CGF_data[0][i];
    d2phi_dlambda2[i] = CGF_data[3][i];
  }
  
  int n_delta = constants::N_delta_values_for_PDFs;
  double var = interpolate_neville_aitken(0.0, &lambda_values, &d2phi_dlambda2, constants::order_of_interpolation);
  double skew = interpolate_neville_aitken_derivative(0.0, &lambda_values, &d2phi_dlambda2, constants::order_of_interpolation);
  double lognormal_shift = lognormal_tools::get_delta0(var, skew);
  double var_Gauss = log(1.0+var/lognormal_shift/lognormal_shift);
  double mean_Gauss = -0.5*var_Gauss;
  double std_Gauss = sqrt(var_Gauss);
  double delta_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double delta_max = lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0);
  double ddelta = (delta_max-delta_min)/double(n_delta);
  
  vector<double> delta_values(n_delta, 0.0);
  for(int i = 0; i < n_delta; i++){
    delta_values[i] = delta_min + double(i)*ddelta;
  }
  
  return this->compute_PDF_from_CGF(&delta_values, &CGF_data);
  
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
      n_lambda = i;
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
  vector<double> phi_prime_prime_prime_values(n_lambda, 0.0);
  
  for(int i = 0; i < n_lambda; i++){
    lambda_values[i] = phi_data[0][i];
    phi_values[i] = phi_data[1][i];
    phi_prime_values[i] = phi_data[2][i];
    phi_prime_prime_values[i] = phi_data[3][i];
    phi_prime_prime_prime_values[i] = phi_data[4][i];
  }
  
  int n_delta = constants::N_delta_values_for_PDFs;
  double var = interpolate_neville_aitken(0.0, &lambda_values, &phi_prime_prime_values, constants::order_of_interpolation);
  double skew = interpolate_neville_aitken_derivative(0.0, &lambda_values, &phi_prime_prime_values, constants::order_of_interpolation);
  double lognormal_shift = lognormal_tools::get_delta0(var, skew);
  double var_Gauss = log(1.0+var/lognormal_shift/lognormal_shift);
  double mean_Gauss = -0.5*var_Gauss;
  double std_Gauss = sqrt(var_Gauss);
  double delta_min = max(-1.0, lognormal_shift*(exp(mean_Gauss-5.0*std_Gauss)-1.0));
  double delta_max = min(lognormal_shift*(exp(mean_Gauss+5.0*std_Gauss)-1.0), delta_c);
  
  /*
   * ISSUE: still avoiding delta=delta_max here, but just for conformity since this was 
   * necessary in compute_LOS_projected_PDF.
   */
  double ddelta = (delta_max-delta_min)/double(n_delta);
  double delta, lambda, phi, phi_prime_prime, phi_prime_prime_prime, phi_4;
  
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
  double next_to_leading_correction;
  cout << "Computing 2D PDF:\n";
  for(int i = 0; i < n_delta; i++){
    delta = delta_min + double(i)*ddelta;
    lambda = interpolate_neville_aitken(delta, &phi_prime_values, &lambda_values, constants::order_of_interpolation);
    phi = interpolate_neville_aitken(delta, &phi_prime_values, &phi_values, constants::order_of_interpolation);
    phi_prime_prime = interpolate_neville_aitken(delta, &phi_prime_values, &phi_prime_prime_values, constants::order_of_interpolation);
    phi_prime_prime_prime = interpolate_neville_aitken(lambda, &lambda_values, &phi_prime_prime_prime_values, constants::order_of_interpolation);
    phi_4 = interpolate_neville_aitken_derivative(lambda, &lambda_values, &phi_prime_prime_prime_values, constants::order_of_interpolation);
    PDF_data[0][i] = delta;
    PDF_data[1][i] = exp(-lambda*delta + phi)/sqrt(constants::pi2*phi_prime_prime);
    
    //next_to_leading_correction = 1.0 + (3.0*phi_4*phi_prime_prime - 5.0*pow(phi_prime_prime_prime, 2))/(24.0*pow(phi_prime_prime, 3));
    next_to_leading_correction = 1.0 + 1.0*phi_4/(8.0*pow(phi_prime_prime, 2));
    next_to_leading_correction -= 5.0*pow(phi_prime_prime_prime, 2)/(24.0*pow(phi_prime_prime, 3));
    if(next_to_leading_correction < 1.1 && next_to_leading_correction > 0.9)
      PDF_data[1][i] *= next_to_leading_correction;
    //else
    //  PDF_data[1][i] = 0.0;
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
  
  vector<double> delta_values_temp(0, 0.0);
  vector<double> F_values(0, 0.0);
  vector<double> F_prime_values(0, 0.0);
  vector<double> F_prime_prime_values(0, 0.0);

  this->return_delta_NL_of_delta_L_and_dF_ddelta_2D(e, &delta_values_temp, &F_values, &F_prime_values, &F_prime_prime_values);
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
  double e;
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
  
  double D_11;
  double D_22;
  double mu;
  double one_plus_mu;
  double vNL, vL, dlnvL_dlnR, S_3;
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time_of_nonzero_kernel; t++){
    int i = indeces_of_nonzero_kernel[t];
    w = w_values_bin_center[i];
    e = eta_0-w;
    R = w*theta;
    
    CGF_data = this->compute_phi_tilde_of_lambda_2D(e, R, f_NL, var_NL_rescale);
    tau_values = CGF_data[0];
    lambda_values = CGF_data[2];
    phi_values = CGF_data[3];
    phi_prime_values = CGF_data[1];
    for(int d = 0; d < tau_values.size(); d++){
      tau_values[d] /= sqrt(CGF_data[4][d]);
    }
    this->compute_polynomial_coefficients_from_CGF(e, R, f_NL, var_NL_rescale, lambda_values, tau_values, phi_values, phi_prime_values, &coefficients_phi_of_lambda[t], &coefficients_phi_prime_of_lambda[t]);
    
    
    D_11 = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_of_delta, constants::order_of_interpolation);
    D_22 = interpolate_neville_aitken(e, &this->eta, &this->Newtonian_growth_factor_second_order, constants::order_of_interpolation);
    mu = 1.0 - D_22/D_11/D_11;
    one_plus_mu = (1.0+mu);
    this->current_P_L = this->P_L(e);
    this->current_P_NL = this->P_NL(e);
    vNL = variance_of_matter_within_R_NL_2D(R)*var_NL_rescale;
    vL = this->variance_of_matter_within_R_2D(R);
    dlnvL_dlnR = this->dvariance_of_matter_within_R_dR_2D(R)/vL*R;
    
    S_3 = 3.0*one_plus_mu + 1.5*dlnvL_dlnR;
    
    coefficients_phi_prime_of_lambda[t][0] = 0.0;
    coefficients_phi_prime_of_lambda[t][1] = vNL;
    coefficients_phi_prime_of_lambda[t][2] = 0.5*S_3*vNL*vNL;
    
    coefficients_phi_of_lambda[t][0] = 0.0;
    coefficients_phi_of_lambda[t][1] = 0.0;
    coefficients_phi_of_lambda[t][2] = 0.5*coefficients_phi_prime_of_lambda[t][1];
    coefficients_phi_of_lambda[t][3] = coefficients_phi_prime_of_lambda[t][2]/3.0;
    
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
 * return_LOS_integrated_phi_of_lambda_incl_CMB_kappa
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::return_LOS_integrated_phi_of_lambda_incl_CMB_kappa(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values, vector<vector<double> > *phi_data_delta, vector<vector<double> > *phi_data_kappa, vector<vector<double> > *phi_grid, vector<vector<double> > *dphi_dldelta_grid, vector<vector<double> > *dphi_dlkappa_grid, vector<vector<double> > *d2phi_dldelta2_grid, vector<vector<double> > *d3phi_dldelta3_grid, vector<vector<double> > *d3phi_dldelta2_dlkappa_grid, vector<vector<double> > *d3phi_dldelta_dlkappa2_grid, vector<vector<double> > *d4phi_dldelta4_grid, vector<vector<double> > *d4phi_dldelta2_dlkappa2_grid, vector<vector<double> > *d2phi_dldelta_dlkappa_grid, vector<vector<double> > *d2phi_dlkappa2_grid, vector<vector<int> > *grid_mask){
  
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
  z = 1.0/this->a_at_eta(eta_0 - w)-1.0;
  while(w < w_last_scattering){ // 2.35 is for Buzzard
    w_values_extended_to_last_scattering.push_back(w);
    kernel_values_extended_to_last_scattering.push_back(0.0);
    n_time = w_values_extended_to_last_scattering.size()-1;
    dw = min(constants::maximal_dw, w_last_scattering - w);
    w = w + dw;
    z = 1.0/this->a_at_eta(eta_0 - w)-1.0;
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
  double z_min, z_max;
  double norm = 0.0;
  for(int i = 0; i < n_time; i++){
    w = w_values_bin_center[i];
    w_min = w_values_extended_to_last_scattering[i];
    w_max = w_values_extended_to_last_scattering[i+1];
    z_min = 1.0/this->a_at_eta(eta_0-w_min)-1.0;
    z_max = 1.0/this->a_at_eta(eta_0-w_max)-1.0;
    
    a = this->a_at_eta(eta_0-w);
    z = 1.0/a-1.0;
    CMB_lensing_kernel[i] = 1.5*this->return_Omega_m()*w*(w_last_scattering-w)/w_last_scattering/a;
  }
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  vector<vector<double> > y_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_prime_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_prime_prime_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_3_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > phi_4_values(n_time, vector<double>(n_lambda, 0.0));
  vector<vector<double> > dummy_data;
  
  double y_min = 0.0;
  double y_max = 0.0;
  vector<double>::iterator y_min_iterator;
  vector<double>::iterator y_max_iterator;
  int time_index_y_max, lambda_index_y_max, lambda_index_y_min;
    
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time; t++){
    w = w_values_bin_center[t];
    eta = eta_0-w;
    R = w*theta;
    dummy_data = compute_phi_tilde_of_lambda_2D(eta, R, f_NL, var_NL_rescale);
    y_values[t] = dummy_data[0];
    phi_values[t] = dummy_data[1];
    phi_prime_values[t] = dummy_data[2];
    phi_prime_prime_values[t] = dummy_data[3];
    phi_3_values[t] = dummy_data[4];
    for(int l = 0; l < phi_4_values[t].size(); l++){
      phi_4_values[t][l] = interpolate_neville_aitken_derivative(y_values[t][l], &y_values[t], &phi_3_values[t], constants::order_of_interpolation);
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
    phi_3_values[t] = vector<double>(&phi_3_values[t][0], &phi_3_values[t][lambda_index_y_max]+1);
    phi_4_values[t] = vector<double>(&phi_4_values[t][0], &phi_4_values[t][lambda_index_y_max]+1);
    
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
  
  int N_ldelta = (*phi_data_delta)[0].size();
  int N_lkappa = (*phi_data_kappa)[0].size();
  (*phi_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*dphi_dldelta_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*dphi_dlkappa_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*d2phi_dldelta2_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*d3phi_dldelta3_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*d3phi_dldelta2_dlkappa_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*d3phi_dldelta_dlkappa2_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*d4phi_dldelta4_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
  (*d4phi_dldelta2_dlkappa2_grid) = vector<vector<double> >(N_ldelta, vector<double>(N_lkappa, 0.0));
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
  double d3phi_dldelta3;
  double d3phi_dldelta2_dlkappa;
  double d3phi_dldelta_dlkappa2;
  double d4phi_dldelta4;
  double d4phi_dldelta2_dlkappa2;
  double d2phi_dldelta_dlkappa;
  double d2phi_dlkappa2;
  
  double phi_prime_integrand;
  double phi_prime_prime_integrand;
  vector<double> * pointer_to_y_values;
  
  cout << "Projecting 2D CGF:\n";

  
  for(int d = 0; d < N_ldelta; d++){
    cout << "\r" << int(10000.0*double(d)/double(N_ldelta))/100.0 << "%";
    cout.flush();
    ldelta = (*phi_data_delta)[0][d];
    for(int k = 0; k < N_lkappa; k++){
      lkappa = (*phi_data_kappa)[0][k];
      
      
      phi = interpolate_neville_aitken(lkappa, &phi_data_kappa_nonoverlap[0], &phi_data_kappa_nonoverlap[1], constants::order_of_interpolation);
      dphi_dldelta = 0.0;
      dphi_dlkappa = interpolate_neville_aitken(lkappa, &phi_data_kappa_nonoverlap[0], &phi_data_kappa_nonoverlap[2], constants::order_of_interpolation);
      d2phi_dldelta2 = 0.0;
      d3phi_dldelta3 = 0.0;
      d3phi_dldelta2_dlkappa = 0.0;
      d3phi_dldelta_dlkappa2 = 0.0;
      d4phi_dldelta4 = 0.0;
      d4phi_dldelta2_dlkappa2 = 0.0;
      d2phi_dldelta_dlkappa = 0.0;
      d2phi_dlkappa2 = interpolate_neville_aitken(lkappa, &phi_data_kappa_nonoverlap[0], &phi_data_kappa_nonoverlap[3], constants::order_of_interpolation);
      
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
        
        phi_prime_integrand = interpolate_neville_aitken(y, pointer_to_y_values, &phi_prime_values[t], constants::order_of_interpolation);
        dphi_dldelta += dw_times_kernel_delta*phi_prime_integrand;
        dphi_dlkappa += dw_times_kernel_kappa*phi_prime_integrand;
        
        lambda_index_y_max = pointer_to_y_values->size() - 1;
        if(y < y_values[t][lambda_index_y_max] && y > y_values[t][0]){
          phi += dw*interpolate_neville_aitken(y, pointer_to_y_values, &phi_values[t], constants::order_of_interpolation);
          phi_prime_prime_integrand = interpolate_neville_aitken(y, pointer_to_y_values, &phi_prime_prime_values[t], constants::order_of_interpolation);
          d2phi_dldelta2 += dw_times_kernel_delta*kernel_delta*phi_prime_prime_integrand;
          d3phi_dldelta3 += dw*pow(kernel_delta, 3)*interpolate_neville_aitken(y, pointer_to_y_values, &phi_3_values[t], constants::order_of_interpolation);
          d3phi_dldelta2_dlkappa += dw_times_kernel_kappa*pow(kernel_delta, 2)*interpolate_neville_aitken(y, pointer_to_y_values, &phi_3_values[t], constants::order_of_interpolation);
          d3phi_dldelta_dlkappa2 += dw_times_kernel_delta*pow(kernel_kappa, 2)*interpolate_neville_aitken(y, pointer_to_y_values, &phi_3_values[t], constants::order_of_interpolation);
          d4phi_dldelta4 += dw*pow(kernel_delta, 4)*interpolate_neville_aitken(y, pointer_to_y_values, &phi_4_values[t], constants::order_of_interpolation);
          d4phi_dldelta2_dlkappa2 += dw*pow(kernel_delta, 2)*pow(kernel_kappa, 2)*interpolate_neville_aitken(y, pointer_to_y_values, &phi_4_values[t], constants::order_of_interpolation);
          d2phi_dldelta_dlkappa += dw_times_kernel_delta*kernel_kappa*phi_prime_prime_integrand;
          d2phi_dlkappa2 += dw_times_kernel_kappa*kernel_kappa*phi_prime_prime_integrand;
          (*grid_mask)[d][k] = 1;
        }
        //else{
        //  (*grid_mask)[d][k] = 0;
        //}
            
      }
          
      (*phi_grid)[d][k] = phi;
      (*dphi_dldelta_grid)[d][k] = dphi_dldelta;
      (*dphi_dlkappa_grid)[d][k] = dphi_dlkappa;
      (*d2phi_dldelta2_grid)[d][k] = d2phi_dldelta2;
      (*d3phi_dldelta3_grid)[d][k] = d3phi_dldelta3;
      (*d3phi_dldelta2_dlkappa_grid)[d][k] = d3phi_dldelta2_dlkappa;
      (*d3phi_dldelta_dlkappa2_grid)[d][k] = d3phi_dldelta_dlkappa2;
      (*d4phi_dldelta4_grid)[d][k] = d4phi_dldelta4;
      (*d4phi_dldelta2_dlkappa2_grid)[d][k] = d4phi_dldelta2_dlkappa2;
      (*d2phi_dldelta_dlkappa_grid)[d][k] = d2phi_dldelta_dlkappa;
      (*d2phi_dlkappa2_grid)[d][k] = d2phi_dlkappa2;
      
    }
  }
  cout << "\r100%         \n";
  
}




/*
 * FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point
 * 
 * 
 * 
 */

void FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point(double theta, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_noise_variance, vector<double> w_values, vector<double> kernel_values, vector<vector<double> > *delta_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid){
  
  
  cout << "Computing projected phi_data:\n";
  cout.flush();
  
  vector<vector<double> > phi_data_delta;
  vector<vector<double> > phi_data_kappa;
  vector<vector<double> > phi_grid;
  vector<vector<double> > dphi_dldelta_grid;
  vector<vector<double> > dphi_dlkappa_grid;
  vector<vector<double> > d2phi_dldelta2_grid;
  vector<vector<double> > d3phi_dldelta3_grid;
  vector<vector<double> > d3phi_dldelta2_dlkappa_grid;
  vector<vector<double> > d3phi_dldelta_dlkappa2_grid;
  vector<vector<double> > d4phi_dldelta4_grid;
  vector<vector<double> > d4phi_dldelta2_dlkappa2_grid;
  vector<vector<double> > d2phi_dldelta_dlkappa_grid;
  vector<vector<double> > d2phi_dlkappa2_grid;
  vector<vector<int> > grid_mask;
  
  this->return_LOS_integrated_phi_of_lambda_incl_CMB_kappa(theta, f_NL, var_NL_rescale, w_values, kernel_values, &phi_data_delta, &phi_data_kappa, &phi_grid, &dphi_dldelta_grid, &dphi_dlkappa_grid, &d2phi_dldelta2_grid, &d3phi_dldelta3_grid, &d3phi_dldelta2_dlkappa_grid, &d3phi_dldelta_dlkappa2_grid, &d4phi_dldelta4_grid, &d4phi_dldelta2_dlkappa2_grid, &d2phi_dldelta_dlkappa_grid, &d2phi_dlkappa2_grid, &grid_mask);
  
  double covariance = interpolate_grid_biquadratic(0.0, 0.0, &phi_data_delta[0], &phi_data_kappa[0], &d2phi_dldelta_dlkappa_grid);
  double dkk = interpolate_grid_biquadratic(0.0, 0.0, &phi_data_delta[0], &phi_data_kappa[0], &d3phi_dldelta_dlkappa2_grid);
  double ddkk = interpolate_grid_biquadratic(0.0, 0.0, &phi_data_delta[0], &phi_data_kappa[0], &d4phi_dldelta2_dlkappa2_grid);
  double ddk = interpolate_grid_biquadratic(0.0, 0.0, &phi_data_delta[0], &phi_data_kappa[0], &d3phi_dldelta2_dlkappa_grid);
  double ddd = interpolate_grid_biquadratic(0.0, 0.0, &phi_data_delta[0], &phi_data_kappa[0], &d3phi_dldelta3_grid);
  double dddd = interpolate_grid_biquadratic(0.0, 0.0, &phi_data_delta[0], &phi_data_kappa[0], &d4phi_dldelta4_grid);
  
  int N_lambda_delta = phi_grid.size();
  int N_lambda_kappa = phi_grid[0].size();
  for(int k = 0; k < N_lambda_kappa; k++){
    phi_data_kappa[1][k] += 0.5*kappa_noise_variance*pow(phi_data_kappa[0][k],2);
    phi_data_kappa[2][k] += kappa_noise_variance*phi_data_kappa[0][k];
    phi_data_kappa[3][k] += kappa_noise_variance;
  }
  
  // on delta-axis we need to compute the PDF on a 5sigma intervall, because
  // shot-noise makes counts-in-cells histogram sensitive to wide range in delta.
  double N_sigma = 5.0;
  double var_delta = interpolate_neville_aitken(0.0, &phi_data_delta[0], &phi_data_delta[3], constants::order_of_interpolation);
  double one_over_var_delta = 1.0/var_delta;
  double std_delta = sqrt(var_delta);
  double skew = interpolate_neville_aitken_derivative(0.0, &phi_data_delta[0], &phi_data_delta[3], constants::order_of_interpolation);
  double lognormal_shift_delta = lognormal_tools::get_delta0(var_delta, skew);
  double var_delta_Gauss = log(1.0+var_delta/lognormal_shift_delta/lognormal_shift_delta);
  double one_over_var_delta_Gauss = 1.0/var_delta_Gauss;
  double mean_delta_Gauss = -0.5*var_delta_Gauss;
  double std_delta_Gauss = sqrt(var_delta_Gauss);
  double delta_min = max(-1.0, lognormal_shift_delta*(exp(mean_delta_Gauss-N_sigma*std_delta_Gauss)-1.0));
  double delta_max = lognormal_shift_delta*(exp(mean_delta_Gauss+N_sigma*std_delta_Gauss)-1.0);
  
  double var_kappa = interpolate_neville_aitken(0.0, &phi_data_kappa[0], &phi_data_kappa[3], constants::order_of_interpolation);
  double one_over_var_kappa = 1.0/var_kappa;
  double std_kappa = sqrt(var_kappa);
  skew = interpolate_neville_aitken_derivative(0.0, &phi_data_kappa[0], &phi_data_kappa[3], constants::order_of_interpolation);
  double lognormal_shift_kappa = lognormal_tools::get_delta0(var_kappa, skew);
  double var_kappa_Gauss = log(1.0+var_kappa/lognormal_shift_kappa/lognormal_shift_kappa);
  double one_over_var_kappa_Gauss = 1.0/var_kappa_Gauss;
  double mean_kappa_Gauss = -0.5*var_kappa_Gauss;
  double std_kappa_Gauss = sqrt(var_kappa_Gauss);
  
  if(kappa_min > lognormal_shift_kappa*(exp(mean_kappa_Gauss-3.0*std_kappa_Gauss)-1.0) || kappa_max < lognormal_shift_kappa*(exp(mean_kappa_Gauss+3.0*std_kappa_Gauss)-1.0)){
    error_handling::general_warning("WARNING in compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point:\nYour intervall [kappa_min, kappa_max] might contain less than 3sigma of probability.");
  }
  
  vector<vector<double> > PDF_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > logPDF_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > delta_Gauss_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > kappa_Gauss_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  
  vector<vector<double> > delta_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > kappa_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > PDF_grid_when_forcing_delta_and_kappa(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  double Lkappa_when_forcing_kappa;
  
  double delta, kappa, delta_Gauss, kappa_Gauss;
  double phi;
  double dphi_dldelta;
  double dphi_dlkappa;
  double d2phi_dldelta2;
  double d3phi_dldelta3;
  double d4phi_dldelta4;
  double d2phi_dldelta_dlkappa;
  double d2phi_dlkappa_dldelta;
  double d2phi_dlkappa2;
  double determinant;
  double next_to_leading_correction;
  
  for(int d = 0; d < N_lambda_delta; d++){
    for(int k = 0; k < N_lambda_kappa; k++){
      // Adding kappa noise variance
      phi_grid[d][k] += 0.5*kappa_noise_variance*pow(phi_data_kappa[0][k],2);
      dphi_dlkappa_grid[d][k] += kappa_noise_variance*phi_data_kappa[0][k];
      d2phi_dlkappa2_grid[d][k] += kappa_noise_variance;
      
      phi = phi_grid[d][k];
      dphi_dldelta = dphi_dldelta_grid[d][k];
      dphi_dlkappa = dphi_dlkappa_grid[d][k];
      d2phi_dldelta2 = d2phi_dldelta2_grid[d][k];
      d3phi_dldelta3 = d3phi_dldelta3_grid[d][k];
      d4phi_dldelta4 = d4phi_dldelta4_grid[d][k];
      d2phi_dldelta_dlkappa = d2phi_dldelta_dlkappa_grid[d][k];
      d2phi_dlkappa_dldelta = d2phi_dldelta_dlkappa;
      d2phi_dlkappa2 = d2phi_dlkappa2_grid[d][k];
      determinant = d2phi_dldelta2*d2phi_dlkappa2;
      determinant -= d2phi_dlkappa_dldelta*d2phi_dldelta_dlkappa;
      
      if(determinant > 0.0 && d2phi_dldelta2 > 0.0 && d2phi_dlkappa2 > 0.0 && grid_mask[d][k] == 1){
        PDF_on_lambda_grid[d][k] = exp(-dphi_dldelta*phi_data_delta[0][d]-dphi_dlkappa*phi_data_kappa[0][k]+phi);
        PDF_on_lambda_grid[d][k] /= constants::pi2*sqrt(determinant);
        next_to_leading_correction = 1.0 + d4phi_dldelta4/(8.0*pow(d2phi_dldelta2, 2));
        next_to_leading_correction -= 5.0*pow(d3phi_dldelta3, 2)/(24.0*pow(d2phi_dldelta2, 3));
        //if(next_to_leading_correction < 1.3 && next_to_leading_correction > 0.7){
        if(next_to_leading_correction < 1.1 && next_to_leading_correction > 0.9){
          PDF_on_lambda_grid[d][k] *= next_to_leading_correction;
        }
        //else{
          //PDF_on_lambda_grid[d][k] = 0.0;
          //grid_mask[d][k] = 0;
        //}
      }
      if(determinant <= 0.0 || d2phi_dldelta2 <= 0.0 || d2phi_dlkappa2 <= 0.0){
        grid_mask[d][k] = 0;
      }
      
      //if(dphi_dldelta>0.0 && dphi_dlkappa > 0.0 && PDF_on_lambda_grid[d-1][k] <= 0.0){
      if(dphi_dldelta>0.0 && dphi_dlkappa > 0.0 && PDF_on_lambda_grid[d][k] <= 0.0){
        PDF_on_lambda_grid[d][k] = 0.0;
      }
      //if(dphi_dldelta>-std_delta && dphi_dlkappa > -std_kappa && PDF_on_lambda_grid[d][k-1] <= 0.0){
      //  PDF_on_lambda_grid[d][k] = 0.0;
      //}
      
    }
  }
  
  vector<vector<double> > dphi_dlkappa_cutout(0);
  vector<vector<double> > dphi_dldelta_cutout(0);
  vector<vector<double> > PDF_cutout(0);
  
  
  for(int d = 0; d < N_lambda_delta; d++){
    PDF_cutout.push_back(vector<double>(0));
    dphi_dlkappa_cutout.push_back(vector<double>(0));
    dphi_dldelta_cutout.push_back(vector<double>(0));
    for(int k = 0; k < N_lambda_kappa; k++){
      if(grid_mask[d][k] > 0){
        PDF_cutout[d].push_back(PDF_on_lambda_grid[d][k]);
        dphi_dlkappa_cutout[d].push_back(dphi_dlkappa_grid[d][k]);
        dphi_dldelta_cutout[d].push_back(dphi_dldelta_grid[d][k]);
      }
    }
  }
  
  
  vector<vector<double> > cutout_from_delta_grid_when_forcing_kappa(N_lambda_kappa, vector<double>(0));
  vector<vector<double> > cutout_from_PDF_grid_when_forcing_kappa(N_lambda_kappa, vector<double>(0));
  
  for(int d = 0; d < N_lambda_delta; d++){
    for(int k = 0; k < N_lambda_kappa; k++){
      delta_on_lambda_grid[d][k] = phi_data_delta[2][d];
      kappa_on_lambda_grid[d][k] = phi_data_kappa[2][k];
      if(dphi_dlkappa_cutout[d].size()>0){
        if(kappa_on_lambda_grid[d][k] >= dphi_dlkappa_cutout[d][0] && kappa_on_lambda_grid[d][k] <= dphi_dlkappa_cutout[d][dphi_dlkappa_cutout[d].size()-1]){
          cutout_from_delta_grid_when_forcing_kappa[k].push_back(interpolate_neville_aitken(kappa_on_lambda_grid[d][k], &dphi_dlkappa_cutout[d], &dphi_dldelta_cutout[d], constants::order_of_interpolation));
          cutout_from_PDF_grid_when_forcing_kappa[k].push_back(interpolate_neville_aitken(kappa_on_lambda_grid[d][k], &dphi_dlkappa_cutout[d], &PDF_cutout[d], constants::order_of_interpolation));
        }
      }
    }
  }
  
  for(int d = 0; d < N_lambda_delta; d++){
    for(int k = 0; k < N_lambda_kappa; k++){
      if(cutout_from_delta_grid_when_forcing_kappa[k].size()>0){
        if(delta_on_lambda_grid[d][k] >= cutout_from_delta_grid_when_forcing_kappa[k][0] && delta_on_lambda_grid[d][k] <= cutout_from_delta_grid_when_forcing_kappa[k][cutout_from_delta_grid_when_forcing_kappa[k].size()-1]){
          PDF_grid_when_forcing_delta_and_kappa[d][k] = max(0.0,interpolate_neville_aitken(delta_on_lambda_grid[d][k], &cutout_from_delta_grid_when_forcing_kappa[k], &cutout_from_PDF_grid_when_forcing_kappa[k], constants::order_of_interpolation));
        }
      }
      
      
      if(delta_on_lambda_grid[d][k]>0.0 && kappa_on_lambda_grid[d][k] > 0.0 && PDF_grid_when_forcing_delta_and_kappa[d-1][k] <= 0.0){
        PDF_grid_when_forcing_delta_and_kappa[d][k] = 0.0;
      }
      //if(delta_on_lambda_grid[d][k]>-std_delta && kappa_on_lambda_grid[d][k] > -std_kappa && PDF_grid_when_forcing_delta_and_kappa[d][k-1] <= 0.0){
      //  PDF_grid_when_forcing_delta_and_kappa[d][k] = 0.0;
      //}
      
    }
  }
  
  
  vector<vector<int> > mask_for_biquadratic_interpolation = grid_mask;
  
  
  for(int d = 0; d < N_lambda_delta; d++){
    for(int k = 0; k < N_lambda_kappa; k++){
      if(PDF_grid_when_forcing_delta_and_kappa[d][k] > 0.0){
        mask_for_biquadratic_interpolation[d][k] = 1;
      }
      else{
        mask_for_biquadratic_interpolation[d][k] = 0;
      }
    }
  }
  
  int n_delta = constants::N_delta_values_for_PDFs;
  int n_kappa = constants::N_delta_values_for_PDFs;
  
  (*delta_grid) = vector<vector<double> >(n_delta, vector<double>(n_kappa, 0.0));
  (*kappa_grid) = vector<vector<double> >(n_delta, vector<double>(n_kappa, 0.0));
  (*PDF_grid) = vector<vector<double> >(n_delta, vector<double>(n_kappa, 0.0));
  
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
  
  FILE *F = fopen("Data/test_PDF_on_physical_grid.dat", "w");
  fclose(F);
  fstream out;
  out.open("Data/test_PDF_on_physical_grid.dat");
  out << scientific << setprecision(10);
  
  

  for(int d = 0; d < n_delta; d++){
    delta = delta_values[d];
    delta_Gauss = log(1.0 + delta/lognormal_shift_delta);
    for(int k = 0; k < n_kappa; k++){
      kappa = kappa_values[k];
      kappa_Gauss = log(1.0 + kappa/lognormal_shift_kappa);
      
      (*delta_grid)[d][k] = delta;
      (*kappa_grid)[d][k] = kappa;
      // bi-quadratic interpolation can become unstable (and is unnecessary) in the tails
      if(delta_Gauss/std_delta_Gauss < 2.0 && kappa_Gauss/std_kappa_Gauss < 2.0)
        (*PDF_grid)[d][k] = interpolate_grid_biquadratic(delta, kappa, &phi_data_delta[2], &phi_data_kappa[2], &PDF_grid_when_forcing_delta_and_kappa);
      else
        (*PDF_grid)[d][k] = interpolate_grid_bilinear(delta, kappa, &phi_data_delta[2], &phi_data_kappa[2], &PDF_grid_when_forcing_delta_and_kappa);
      
      if((*PDF_grid)[d][k] < 0.0) (*PDF_grid)[d][k] = 0.0;
      if(delta < phi_data_delta[2][0] || delta > phi_data_delta[2][N_lambda_delta-1]) (*PDF_grid)[d][k] = 0.0;
      if(kappa < phi_data_kappa[2][0] || kappa > phi_data_kappa[2][N_lambda_kappa-1]) (*PDF_grid)[d][k] = 0.0;
      
      
      //if(delta>0.0 && kappa > 0.0 && (*PDF_grid)[d-1][k] <= 0.0){
      //  (*PDF_grid)[d][k] = 0.0;
      //}
      //if(delta>-std_delta && kappa > -std_kappa && (*PDF_grid)[d][k-1] <= 0.0){
      //  (*PDF_grid)[d][k] = 0.0;
      //}
      
      out << (*PDF_grid)[d][k] << setw(20);
      
    }
    out << '\n';
  }
  out.close();
  
  double norm = 0.0;
  double maxi = 0.0;
  int dmaxi = 0;
  int kmaxi = 0;
  for(int d = 0; d < n_delta; d++){
    for(int k = 0; k < n_kappa; k++){
      norm += (*PDF_grid)[d][k];
      if((*PDF_grid)[d][k] > maxi){
        maxi = (*PDF_grid)[d][k];
        dmaxi = d;
        kmaxi = k;
      }
    }
  }
  norm *= ddelta*dkappa;
  
  cout << scientific << setprecision(5);
  cout << "kappa variance = " << var_kappa  << "\n";
  cout << "delta variance = " << var_delta  << "\n";
  cout << "covariance = " << covariance  << "\n";
  cout << "dkk = " << dkk  << "\n";
  cout << "ddkk = " << ddkk  << "\n";
  cout << "ddd = " << ddd  << "\n";
  cout << "dddd = " << dddd  << "\n";
  cout << "ddk = " << ddk  << "\n";
  cout << "# norm:\n";
  cout << "# " << norm << '\n';
  cout << "# max:\n";
  cout << "# " << dmaxi << "  " << kmaxi << "  " << maxi << '\n';
  cout << "# P[250,250]:\n";
  cout << "# " << (*PDF_grid)[250][250] << '\n';
  cout << "# delta bounds:\n";
  cout << "# " << delta_min << "  " << delta_max << '\n';
  cout << "# " << phi_data_delta[2][0] << "  " << phi_data_delta[2][N_lambda_delta-1] << '\n';
  cout << "# kappa bounds:\n";
  cout << "# " << kappa_min << "  " << kappa_max << '\n';
  cout << "# " << phi_data_kappa[2][0] << "  " << phi_data_kappa[2][N_lambda_kappa-1] << '\n';
  
}



/*******************************************************************************************************************************************************
 * return_LOS_integrated_phi_of_lambda_incl_kappa
 * 
 * 
*******************************************************************************************************************************************************/

void FlatInhomogeneousUniverseLCDM::return_LOS_integrated_phi_of_lambda_incl_kappa(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<double> kernel_values, vector<double> lensing_kernel_values, vector<vector<double> > *phi_data_delta, vector<vector<double> > *phi_data_kappa, vector<vector<double> > *phi_grid, vector<vector<double> > *dphi_dldelta_grid, vector<vector<double> > *dphi_dlkappa_grid, vector<vector<double> > *d2phi_dldelta2_grid, vector<vector<double> > *d2phi_dldelta_dlkappa_grid, vector<vector<double> > *d2phi_dlkappa2_grid, vector<vector<int> > *grid_mask){
  
  int n_lambda = this->delta_values_for_cylindrical_collapse.size();
  int n_time = w_values.size()-1;
  double z, a, eta, eta_0, w, dw, w_last_scattering, R;
  eta_0 = this->eta_at_a(1.0);
  
  vector<double> w_values_bin_center(n_time, 0.0);
  for(int i = 0; i < n_time; i++){
    w_values_bin_center[i] = 0.5*(w_values[i+1]+w_values[i]);
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
  
  cout << "computing CGF grid & cutting out 1st branch\n";
  
  for(int t = 0; t < n_time; t++){
    w = w_values_bin_center[t];
    eta = eta_0-w;
    R = w*theta;
    dummy_data = compute_phi_tilde_of_lambda_2D(eta, R, f_NL, var_NL_rescale);
    y_values[t] = dummy_data[0];
    phi_values[t] = dummy_data[1];
    phi_prime_values[t] = dummy_data[2];
    phi_prime_prime_values[t] = dummy_data[3];
    
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
  
  cout << "computing ranges of lambda_delta and lambda_kappa & resizing phi_grids accordingly\n";
  
  vector<int> indeces_of_nonzero_kernel(0,0);
  vector<double> lensing_kernel_nonoverlap = lensing_kernel_values;
  for(int t = 0; t < n_time; t++){
    if(kernel_values[t] > 0.0){
      indeces_of_nonzero_kernel.push_back(t);
      lensing_kernel_nonoverlap[t] = 0.0;
    }
  }
  
  (*phi_data_delta) = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, var_NL_rescale, w_values, kernel_values);
  (*phi_data_kappa) = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, var_NL_rescale, w_values, lensing_kernel_values);
  vector<vector<double> > phi_data_kappa_nonoverlap = this->return_LOS_integrated_phi_of_lambda(theta, f_NL, var_NL_rescale, w_values, lensing_kernel_nonoverlap);
  
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
  
  cout << "Projecting 2D CGF:\n";
  
  for(int d = 0; d < N_ldelta; d++){
    cout << "\r" << int(10000.0*double(d)/double(N_ldelta))/100.0 << "%";
    cout.flush();
    
    ldelta = (*phi_data_delta)[0][d];
    for(int k = 0; k < N_lkappa; k++){
      lkappa = (*phi_data_kappa)[0][k];
      
      
      phi = interpolate_neville_aitken(lkappa, &phi_data_kappa_nonoverlap[0], &phi_data_kappa_nonoverlap[1], constants::order_of_interpolation);
      dphi_dldelta = 0.0;
      dphi_dlkappa = interpolate_neville_aitken(lkappa, &phi_data_kappa_nonoverlap[0], &phi_data_kappa_nonoverlap[2], constants::order_of_interpolation);
      d2phi_dldelta2 = 0.0;
      d2phi_dldelta_dlkappa = 0.0;
      d2phi_dlkappa2 = interpolate_neville_aitken(lkappa, &phi_data_kappa_nonoverlap[0], &phi_data_kappa_nonoverlap[3], constants::order_of_interpolation);
      
      for(int i = 0; i < N_nonzero_kernel; i++){
        t = indeces_of_nonzero_kernel[i];
        dw = w_values[t+1]-w_values[t];
        kernel_delta = kernel_values[t];
        dw_times_kernel_delta = kernel_delta*dw;
        kernel_kappa = lensing_kernel_values[t];
        dw_times_kernel_kappa = kernel_kappa*dw;
        y = kernel_delta*ldelta;
        y += kernel_kappa*lkappa;
        
        pointer_to_y_values = &y_values[t];
        
        phi_prime_integrand = interpolate_neville_aitken(y, pointer_to_y_values, &phi_prime_values[t], constants::order_of_interpolation);
        dphi_dldelta += dw_times_kernel_delta*phi_prime_integrand;
        dphi_dlkappa += dw_times_kernel_kappa*phi_prime_integrand;
        
        lambda_index_y_max = pointer_to_y_values->size() - 1;
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
      (*dphi_dldelta_grid)[d][k] = dphi_dldelta;
      (*dphi_dlkappa_grid)[d][k] = dphi_dlkappa;
      (*d2phi_dldelta2_grid)[d][k] = d2phi_dldelta2;
      (*d2phi_dldelta_dlkappa_grid)[d][k] = d2phi_dldelta_dlkappa;
      (*d2phi_dlkappa2_grid)[d][k] = d2phi_dlkappa2;
      
    }
  }
  cout << "\r100%         \n";
  
}




/*
 * FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF_incl_kappa_saddle_point
 * 
 * 
 * 
 */



void FlatInhomogeneousUniverseLCDM::compute_LOS_projected_PDF_incl_kappa_saddle_point(double theta, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, vector<double> w_values, vector<double> kernel_values, vector<double> lensing_kernel_values, vector<vector<double> > *delta_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid){
  
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
  
  this->return_LOS_integrated_phi_of_lambda_incl_kappa(theta, f_NL, var_NL_rescale, w_values, kernel_values, lensing_kernel_values, &phi_data_delta, &phi_data_kappa, &phi_grid, &dphi_dldelta_grid, &dphi_dlkappa_grid, &d2phi_dldelta2_grid, &d2phi_dldelta_dlkappa_grid, &d2phi_dlkappa2_grid, &grid_mask);
  
  cout << "kappa variance = " << interpolate_neville_aitken(0.0, &phi_data_kappa[0], &phi_data_kappa[3], constants::order_of_interpolation)  << "\n";
  
  int N_lambda_delta = phi_grid.size();
  int N_lambda_kappa = phi_grid[0].size();
  
  vector<vector<double> > PDF_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > logPDF_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > delta_Gauss_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > kappa_Gauss_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  
  vector<vector<double> > delta_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > kappa_on_lambda_grid(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  vector<vector<double> > delta_grid_when_forcing_kappa_transposed(N_lambda_kappa, vector<double>(N_lambda_delta,0.0));
  vector<vector<double> > PDF_grid_when_forcing_kappa_transposed(N_lambda_kappa, vector<double>(N_lambda_delta,0.0));
  vector<vector<double> > PDF_grid_when_forcing_delta_and_kappa(N_lambda_delta, vector<double>(N_lambda_kappa,0.0));
  double Lkappa_when_forcing_kappa;
  
  double delta, kappa, delta_Gauss, kappa_Gauss;
  double phi;
  double dphi_dldelta;
  double dphi_dlkappa;
  double d2phi_dldelta2;
  double d2phi_dldelta_dlkappa;
  double d2phi_dlkappa_dldelta;
  double d2phi_dlkappa2;
  double determinant;
  for(int d = 0; d < N_lambda_delta; d++){
    for(int k = 0; k < N_lambda_kappa; k++){
      phi = phi_grid[d][k];
      dphi_dldelta = dphi_dldelta_grid[d][k];
      dphi_dlkappa = dphi_dlkappa_grid[d][k];
      d2phi_dldelta2 = d2phi_dldelta2_grid[d][k];
      d2phi_dldelta_dlkappa = d2phi_dldelta_dlkappa_grid[d][k];
      d2phi_dlkappa_dldelta = d2phi_dldelta_dlkappa;
      d2phi_dlkappa2 = d2phi_dlkappa2_grid[d][k];
      determinant = d2phi_dldelta2*d2phi_dlkappa2;
      determinant -= d2phi_dlkappa_dldelta*d2phi_dldelta_dlkappa;
      
      if(determinant > 0.0 && grid_mask[d][k] == 1){
        PDF_on_lambda_grid[d][k] = exp(-dphi_dldelta*phi_data_delta[0][d]-dphi_dlkappa*phi_data_kappa[0][k]+phi);
        PDF_on_lambda_grid[d][k] /= constants::pi2*sqrt(determinant);
      }
      if(determinant <= 0.0){
        grid_mask[d][k] = 0;
      }
    }
  }
  
  
  for(int d = 0; d < N_lambda_delta; d++){
    for(int k = 0; k < N_lambda_kappa; k++){
      delta_on_lambda_grid[d][k] = phi_data_delta[2][d];
      kappa_on_lambda_grid[d][k] = phi_data_kappa[2][k];
      Lkappa_when_forcing_kappa = interpolate_neville_aitken(kappa_on_lambda_grid[d][k], &dphi_dlkappa_grid[d], &phi_data_kappa[0], constants::order_of_interpolation);
      delta_grid_when_forcing_kappa_transposed[k][d] = interpolate_neville_aitken(Lkappa_when_forcing_kappa, &phi_data_kappa[0], &dphi_dldelta_grid[d], constants::order_of_interpolation);
      PDF_grid_when_forcing_kappa_transposed[k][d] = interpolate_neville_aitken(Lkappa_when_forcing_kappa, &phi_data_kappa[0], &PDF_on_lambda_grid[d], constants::order_of_interpolation);
    }
  }
  
  
  for(int d = 0; d < N_lambda_delta; d++){
    for(int k = 0; k < N_lambda_kappa; k++){
      // d_values = delta_grid_when_forcing_kappa_transposed[k][:]
      // p_values = PDF_grid_when_forcing_kappa_transposed[k][:]
      if(delta_on_lambda_grid[d][k] >= delta_grid_when_forcing_kappa_transposed[k][0] && delta_on_lambda_grid[d][k] <= delta_grid_when_forcing_kappa_transposed[k][N_lambda_delta-1])
        PDF_grid_when_forcing_delta_and_kappa[d][k] = max(0.0,interpolate_neville_aitken(delta_on_lambda_grid[d][k], &delta_grid_when_forcing_kappa_transposed[k], &PDF_grid_when_forcing_kappa_transposed[k], constants::order_of_interpolation));
    }
  }
  
  // on delta-axis we need to compute the PDF on a 5sigma intervall, because
  // shot-noise makes counts-in-cells histogram sensitive to wide range in delta.
  double N_sigma = 5.0;
  double var_delta = interpolate_neville_aitken(0.0, &phi_data_delta[0], &phi_data_delta[3], constants::order_of_interpolation);
  double one_over_var_delta = 1.0/var_delta;
  double std_delta = sqrt(var_delta);
  double skew = interpolate_neville_aitken_derivative(0.0, &phi_data_delta[0], &phi_data_delta[3], constants::order_of_interpolation);
  double lognormal_shift_delta = lognormal_tools::get_delta0(var_delta, skew);
  double var_delta_Gauss = log(1.0+var_delta/lognormal_shift_delta/lognormal_shift_delta);
  double one_over_var_delta_Gauss = 1.0/var_delta_Gauss;
  double mean_delta_Gauss = -0.5*var_delta_Gauss;
  double std_delta_Gauss = sqrt(var_delta_Gauss);
  double delta_min = max(-1.0, lognormal_shift_delta*(exp(mean_delta_Gauss-N_sigma*std_delta_Gauss)-1.0));
  double delta_max = lognormal_shift_delta*(exp(mean_delta_Gauss+N_sigma*std_delta_Gauss)-1.0);
  
  double var_kappa = interpolate_neville_aitken(0.0, &phi_data_kappa[0], &phi_data_kappa[3], constants::order_of_interpolation);
  double one_over_var_kappa = 1.0/var_kappa;
  double std_kappa = sqrt(var_kappa);
  skew = interpolate_neville_aitken_derivative(0.0, &phi_data_kappa[0], &phi_data_kappa[3], constants::order_of_interpolation);
  double lognormal_shift_kappa = lognormal_tools::get_delta0(var_kappa, skew);
  double var_kappa_Gauss = log(1.0+var_kappa/lognormal_shift_kappa/lognormal_shift_kappa);
  double one_over_var_kappa_Gauss = 1.0/var_kappa_Gauss;
  double mean_kappa_Gauss = -0.5*var_kappa_Gauss;
  double std_kappa_Gauss = sqrt(var_kappa_Gauss);
  
  if(kappa_min > lognormal_shift_kappa*(exp(mean_kappa_Gauss-3.0*std_kappa_Gauss)-1.0) || kappa_max < lognormal_shift_kappa*(exp(mean_kappa_Gauss+3.0*std_kappa_Gauss)-1.0)){
    error_handling::general_warning("WARNING in compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point:\nYour intervall [kappa_min, kappa_max] might contain less than 3sigma of probability.");
  }
  
  int n_delta = constants::N_delta_values_for_PDFs;
  int n_kappa = constants::N_delta_values_for_PDFs;
  
  (*delta_grid) = vector<vector<double> >(n_delta, vector<double>(n_kappa, 0.0));
  (*kappa_grid) = vector<vector<double> >(n_delta, vector<double>(n_kappa, 0.0));
  (*PDF_grid) = vector<vector<double> >(n_delta, vector<double>(n_kappa, 0.0));
  
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
  
  for(int d = 0; d < n_delta; d++){
    delta = delta_values[d];
    for(int k = 0; k < n_kappa; k++){
      kappa = kappa_values[k];
      
      (*delta_grid)[d][k] = delta;
      (*kappa_grid)[d][k] = kappa;
      (*PDF_grid)[d][k] = interpolate_grid_biquadratic(delta, kappa, &phi_data_delta[2], &phi_data_kappa[2], &PDF_grid_when_forcing_delta_and_kappa);
      if((*PDF_grid)[d][k] < 0.0) (*PDF_grid)[d][k] = 0.0;
      
    }
  }
  
  double norm = 0.0;
  double maxi = 0.0;
  int dmaxi = 0;
  int kmaxi = 0;
  for(int d = 0; d < n_delta; d++){
    for(int k = 0; k < n_kappa; k++){
      norm += (*PDF_grid)[d][k];
      if((*PDF_grid)[d][k] > maxi){
        maxi = (*PDF_grid)[d][k];
        dmaxi = d;
        kmaxi = k;
      }
    }
  }
  norm *= ddelta*dkappa;
  
  cout << scientific << setprecision(5);
  cout << "# norm:\n";
  cout << "# " << norm << '\n';
  cout << "# max:\n";
  cout << "# " << dmaxi << "  " << kmaxi << "  " << maxi << '\n';
  cout << "# P[250,250]:\n";
  cout << "# " << (*PDF_grid)[250][250] << '\n';
  cout << "# delta bounds:\n";
  cout << "# " << delta_min << "  " << delta_max << '\n';
  cout << "# " << phi_data_delta[2][0] << "  " << phi_data_delta[2][N_lambda_delta-1] << '\n';
  cout << "# kappa bounds:\n";
  cout << "# " << kappa_min << "  " << kappa_max << '\n';
  cout << "# " << phi_data_kappa[2][0] << "  " << phi_data_kappa[2][N_lambda_kappa-1] << '\n';
  
}





vector<double> FlatInhomogeneousUniverseLCDM::return_LOS_integrated_C_ells(int l_max, vector<double> w_values, vector<double> kernel_values){
  
  int N_ell = l_max + 1;
  vector<double> C_ells(N_ell, 0.0);
  
  int n_time = w_values.size()-1;
  
  vector<int> indeces_of_nonzero_kernel(0,0);
  
  double eta, eta_0, w, dw, ell, ell_plus_half, T_ell, k;
  eta_0 = this->eta_at_a(1.0);
  
  for(int i = 0; i < n_time; i++){
    if(kernel_values[i] > 0.0){
      indeces_of_nonzero_kernel.push_back(i);
    }
  }
  
  int n_time_of_nonzero_kernel = indeces_of_nonzero_kernel.size();
  
  for(int t = 0; t < n_time_of_nonzero_kernel; t++){
    int i = indeces_of_nonzero_kernel[t];
    w = 0.5*(w_values[i+1]+w_values[i]);
    dw = w_values[i+1]-w_values[i];
    eta = eta_0-w;
    this->current_P_NL = this->P_NL(eta);
    for(int l = 2; l < N_ell; l++){
      ell = float(l);
      ell_plus_half = ell+0.5;
      k = ell_plus_half/w;
      C_ells[l] += pow(kernel_values[i]/w, 2)*interpolate_neville_aitken(log(k), &this->log_wave_numbers, &this->current_P_NL, constants::order_of_interpolation);
    }
  }
  
  for(int l = 2; l < N_ell; l++){
    ell = float(l);
    ell_plus_half = ell+0.5;
    // Limber correction from https://arxiv.org/abs/1611.04954 :
    T_ell = ((ell+2.0)*(ell+1.0)*ell*(ell-1.0))/pow(ell_plus_half, 4.0);
    C_ells[l] *= T_ell;
  }
  
  return C_ells;
  
}





vector<vector<vector<double> > > FlatInhomogeneousUniverseLCDM::return_LOS_integrated_C_ells(int l_max, vector<double> w_values, vector<vector<double> > kernel_values){
  
  int N_ell = l_max + 1;
  int N_fields = kernel_values.size();
  vector<vector<vector<double> > > C_ells(N_fields, vector<vector<double> >(N_fields, vector<double>(N_ell, 0.0)));
  
  int n_time = w_values.size()-1;
  
  vector<int> indeces_of_nonzero_kernel(0,0);
  
  int kernel_positive = 0;
  double eta, eta_0, w, dw, ell, ell_plus_half, T_ell, k;
  eta_0 = this->eta_at_a(1.0);
  
  for(int i = 0; i < n_time; i++){
    kernel_positive = 0;
    for(int k = 0; k < N_fields; k++){
      if(kernel_values[k][i] > 0.0) kernel_positive = 1;
    }
    if(kernel_positive == 1){
      indeces_of_nonzero_kernel.push_back(i);
    }
  }
  
  int n_time_of_nonzero_kernel = indeces_of_nonzero_kernel.size();
  
  cout << "computing C_ells:\n";
  
  for(int t = 0; t < n_time_of_nonzero_kernel; t++){
    int i = indeces_of_nonzero_kernel[t];
    w = 0.5*(w_values[i+1]+w_values[i]);
    dw = w_values[i+1]-w_values[i];
    eta = eta_0-w;
    this->current_P_NL = this->P_NL(eta);
    for(int l = 2; l < N_ell; l++){
      ell = float(l);
      ell_plus_half = ell+0.5;
      k = ell_plus_half/w;
      
      for(int f1 = 0; f1 < N_fields; f1++){
        for(int f2 = f1; f2 < N_fields; f2++){
          C_ells[f1][f2][l] += dw*kernel_values[f1][i]*kernel_values[f2][i]/pow(w, 2)*interpolate_neville_aitken(log(k), &this->log_wave_numbers, &this->current_P_NL, constants::order_of_interpolation);
        }
      }
    }
  }
  
  cout << "Done.\n";
  
  for(int l = 2; l < N_ell; l++){
    ell = float(l);
    ell_plus_half = ell+0.5;
    // Limber correction from https://arxiv.org/abs/1611.04954 :
    T_ell = ((ell+2.0)*(ell+1.0)*ell*(ell-1.0))/pow(ell_plus_half, 4.0);
    for(int f1 = 0; f1 < N_fields; f1++){
      for(int f2 = f1; f2 < N_fields; f2++){
        C_ells[f1][f2][l] *= T_ell;
        C_ells[f2][f1][l] = C_ells[f1][f2][l];
      }
    }
  }
  cout << "Done.\n";
  
  return C_ells;
  
}





vector<vector<vector<double> > > FlatInhomogeneousUniverseLCDM::return_LOS_integrated_3rd_moments(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<vector<double> > kernel_values){
  
  if(f_NL != 0.0){
    cerr << "CAREFUL: f_NL != 0 not yet implemented in return_LOS_integrated_skewness.\n";
  }
  
  int n_time = w_values.size()-1;
  int N_fields = kernel_values.size();
  int kernel_positive = 0;
  double a, e, eta_0, w, dw, R;
  eta_0 = this->eta_at_a(1.0);
  
  vector<int> indeces_of_nonzero_kernel(0,0);
  for(int i = 0; i < n_time; i++){
    kernel_positive = 0;
    for(int k = 0; k < N_fields; k++){
      if(kernel_values[k][i] > 0.0) kernel_positive = 1;
    }
    if(kernel_positive == 1){
      indeces_of_nonzero_kernel.push_back(i);
    }
  }
  
  int n_time_of_nonzero_kernel = indeces_of_nonzero_kernel.size();
  vector<vector<vector<double> > > third_moments(N_fields, vector<vector<double> >(N_fields, vector<double>(N_fields, 0.0)));
  
  double D_11;
  double D_22;
  double mu;
  double one_plus_mu;
  double vNL, vL, dlnvL_dlnR, S_3, skew_cylinder;
  
  cout << "computing 3rd moments:\n";
  
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
    skew_cylinder = S_3*vNL*vNL;
    for(int f1 = 0; f1 < N_fields; f1++){
      for(int f2 = 0; f2 < N_fields; f2++){
        for(int f3 = 0; f3 < N_fields; f3++){
          third_moments[f1][f2][f3] += dw*skew_cylinder*kernel_values[f1][i]*kernel_values[f2][i]*kernel_values[f3][i]; 
        }
      }
    }  
  }
  cout << "Done.\n";
  
  return third_moments;
  
}




vector<vector<double> > FlatInhomogeneousUniverseLCDM::return_LOS_integrated_2nd_moments(double theta, double f_NL, double var_NL_rescale, vector<double> w_values, vector<vector<double> > kernel_values){
  
  if(f_NL != 0.0){
    cerr << "CAREFUL: f_NL != 0 not yet implemented in return_LOS_integrated_skewness.\n";
  }
  
  int n_time = w_values.size()-1;
  int N_fields = kernel_values.size();
  int kernel_positive = 0;
  double a, e, eta_0, w, dw, R;
  eta_0 = this->eta_at_a(1.0);
  
  vector<int> indeces_of_nonzero_kernel(0,0);
  for(int i = 0; i < n_time; i++){
    kernel_positive = 0;
    for(int k = 0; k < N_fields; k++){
      if(kernel_values[k][i] > 0.0) kernel_positive = 1;
    }
    if(kernel_positive == 1){
      indeces_of_nonzero_kernel.push_back(i);
    }
  }
  
  int n_time_of_nonzero_kernel = indeces_of_nonzero_kernel.size();
  vector<vector<double> > second_moments(N_fields, vector<double>(N_fields, 0.0));
  
  double D_11;
  double D_22;
  double mu;
  double one_plus_mu;
  double vNL, vL, dlnvL_dlnR, S_3, skew_cylinder;
  
  cout << "computing 2nd moments:\n";
  
  for(int t = 0; t < n_time_of_nonzero_kernel; t++){
    int i = indeces_of_nonzero_kernel[t];
    w = 0.5*(w_values[i+1]+w_values[i]);
    dw = w_values[i+1]-w_values[i];
    e = eta_0-w;
    R = w*theta;
    this->current_P_NL = this->P_NL(e);
    vNL = variance_of_matter_within_R_NL_2D(R)*var_NL_rescale;
    for(int f1 = 0; f1 < N_fields; f1++){
      for(int f2 = f1; f2 < N_fields; f2++){
        second_moments[f1][f2] += dw*vNL*kernel_values[f1][i]*kernel_values[f2][i];
      }
    }  
  }
  cout << "Done.\n";
  
  for(int f1 = 0; f1 < N_fields; f1++){
    for(int f2 = f1; f2 < N_fields; f2++){
      second_moments[f2][f1] = second_moments[f1][f2];
    }
  } 
  cout << "Done.\n";
  
  return second_moments;
  
}


double FlatInhomogeneousUniverseLCDM::delta_crit(double a, double Delta){
  
  vector<double> delta_L_today_aux;
  vector<double> delta_NL_aux;
  vector<double> delta_NL_prime_aux;
  vector<double> delta_NL_prime_prime_aux;
  
  double e = this->eta_at_a(a);
  double D_growth = this->return_D_of_eta(e);
  this->return_delta_NL_of_delta_L_and_dF_ddelta_3D(e, &delta_L_today_aux, &delta_NL_aux, &delta_NL_prime_aux, &delta_NL_prime_prime_aux);
  double dNL_max = delta_NL_aux[0];
  for(int i = 0; i < delta_NL_aux.size(); i++){
      if(delta_NL_aux[i] > dNL_max) dNL_max = delta_NL_aux[i];
  }
  
  return D_growth*interpolate_neville_aitken(Delta, &delta_NL_aux, &delta_L_today_aux, constants::order_of_interpolation);
  
}
