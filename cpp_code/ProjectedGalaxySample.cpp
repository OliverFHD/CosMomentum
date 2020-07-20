

ProjectedGalaxySample::ProjectedGalaxySample(FlatInhomogeneousUniverseLCDM* universe, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, string n_of_z_input_file) : GalaxySample(universe, b1, b2, a0, a1){
  this->set_n_of_w_data(n_of_z_input_file);
  this->density = density_in_arcmin_squared/pow(constants::pi/(180.0*60.0), 2);
}


ProjectedGalaxySample::~ProjectedGalaxySample(){
}



void ProjectedGalaxySample::set_parameters_projected(double density_in_arcmin_squared, double b1, double b2, double a0, double a1, string n_of_z_input_file){
  this->set_n_of_w_data(n_of_z_input_file);
  this->density = density_in_arcmin_squared/pow(constants::pi/(180.0*60.0), 2);
  this->set_parameters(b1, b2, a0, a1);
}



void ProjectedGalaxySample::set_projected_bias_model_from_br_parametrisation(double b_tilde, double r, double theta_in_arcmin, double f_NL, double var_NL_rescale){
  double A = 2.0*constants::pi*(1.0-cos(theta_in_arcmin*constants::arcmin));
  double N_bar = A*this->density;
  double variance = this->compute_variance_in_angular_tophat(theta_in_arcmin, var_NL_rescale);
  double skewness = this->compute_skewness_in_angular_tophat(theta_in_arcmin, f_NL, var_NL_rescale);
  
  this->set_bias_model_from_br_parametrisation(b_tilde, r, N_bar, variance, skewness);
}


/*
 * ProjectedGalaxySample::set_n_of_w_data
 * 
 * Read in the redshift distribution of the galaxy sample and turn it into a co-moving
 * distance distribution. The textfile "n_of_z_input_file" should have two columns, where
 * column 1 contains the lower-redshift edge of each histogram bin and column 2 contains
 * the redshift histrogram (the histogram doesn't need to be normalised, since normalisation
 * is enforced later on in the code).
 * 
 * ISSUE: there is a potential ambiguity here: the second column of n_of_z_input_file is
 * supposed to be (propotional to) the redshift-density within the redshift bins given by
 * the first column. However, some users might think that the second column should contain
 * the probability within each redshift bin.
 * 
 */
void ProjectedGalaxySample::set_n_of_w_data(string n_of_z_input_file){
  
  vector<vector<double> > nofz_data = read_table_double(n_of_z_input_file);
  
  double eta_0 = this->pointer_to_universe()->eta_at_a(1.0);
  double w_last_scattering = eta_0 - this->pointer_to_universe()->eta_at_a(1.0/(1.0+constants::z_last_scattering));
  int Nz = nofz_data.size();
  double scale_factor;
  double w = 0.0;
  double dw = 0.0;
  double w_input = 0.0;
  double dw_input = 0.0;
  double dz_input = 0.0;
  double n_of_w_input = 0.0;
  
  vector<double> w_dummies(0, 0.0);
  vector<double> n_of_w_dummies(0, 0.0);
  
  // refining the binning of the input file so that dw <= constants::maximal_dw;
  // also, starting the new array at z=0=w, as needed for computation of lensing kernel.
  for(int i = 0; i < Nz; i++){
    dw_input = -w_input; // this line makes sense; do not delete
    scale_factor = 1.0/(1.0+nofz_data[i][0]);
    w_input = eta_0 - this->pointer_to_universe()->eta_at_a(scale_factor);
    if(i > 0){
      dw_input += w_input;
      dz_input = nofz_data[i][0] - nofz_data[i-1][0];
      error_handling::test_value_of_double(dw_input, 0.0, error_handling::LARGER, "nofz_data on n_of_z_input_file not monotonically increasing in z - see ProjectedGalaxySample::set_n_of_w_data(string n_of_z_input_file).");
      n_of_w_input = nofz_data[i-1][1]*dz_input/dw_input;
    }
    while(w < w_input){
      w_dummies.push_back(w);
      if(i == 0)
        n_of_w_dummies.push_back(0.0);
      else
        n_of_w_dummies.push_back(n_of_w_input);
      w += constants::maximal_dw;
    }
    w = w_input;
  }
  w_dummies.push_back(w);
  n_of_w_dummies.push_back(0.0); // n_of_w_dummies[i] is density between w_dummies[i] and w_dummies[i+1], so last value is zero.
  
  Nz = n_of_w_dummies.size();
  double norm = 0.0;
  for(int i = 0; i < Nz-1; i++){
    norm += n_of_w_dummies[i]*(w_dummies[i+1]-w_dummies[i]);
  }
  for(int i = 0; i < Nz; i++){
    n_of_w_dummies[i] /= norm;
  }
  
  
  // copy the new vectors into Class attributes
  // (w_dummies and n_of_w_dummies were needed because vector::push_back can be unefficient in storage, i.e. it can be better to allocate required space in one badge)
  this->w_values = w_dummies;
  this->n_of_w_values = n_of_w_dummies;
  this->lensing_kernel_values = vector<double>(Nz, 0.0);
  
  double scale;
  double w_1, w_2;
  double w_final = this->w_values[Nz-1];
  
  // now calculating the lensing kernel (if ProjectedGalaxySample is used as source galaxy sample)
  for(int i = 0; i < Nz-1; i++){
    w = 0.5*(this->w_values[i]+this->w_values[i+1]);
    w_1 = this->w_values[i+1];
    dw = w_1 - w;
    scale = this->pointer_to_universe()->a_at_eta(eta_0 - w);
    this->lensing_kernel_values[i] += dw*0.5*(w*(w_1-w)/scale/w_1);
    for(int j = i+1; j < Nz-1; j++){
      w_1 = this->w_values[j];
      w_2 = this->w_values[j+1];
      dw = w_2 - w_1;
      //w*(w_final-w)/w_final/a
      this->lensing_kernel_values[i] += dw*0.5*(w*(w_1-w)/w_1);
      this->lensing_kernel_values[i] += dw*0.5*(w*(w_2-w)/w_2);
    }
    this->lensing_kernel_values[i] *= 1.5*this->pointer_to_universe()->return_Omega_m()/scale;
  
    cout << i << "   ";
    cout << this->w_values[i] << "   ";
    cout << this->lensing_kernel_values[i] << "\n";
  }
}


/*
 * ProjectedGalaxySample::compute_variance_in_angular_tophat
 * 
 * Compute variance of the projected density contrast smoothed with a angular top-hat filter of radius theta. The variance allows for a re-scaling wrt the fiducial power spectrum (for var_NL_rescale=1.0 the code uses the Takahashi et al. (2012) version of halofit (see also Smith et al. 2003)).
 * 
 */

double ProjectedGalaxySample::compute_variance_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale){
  return this->pointer_to_universe()->return_LOS_integrated_variance(theta_in_arcmin*constants::arcmin, this->w_values, this->n_of_w_values, var_NL_rescale);
}

/*
 * ProjectedGalaxySample::compute_skewness_in_angular_tophat
 * 
 * Compute 3rd central moment of the projected density contrast smoothed with a angular top-hat filter of radius theta. The variance allows for a re-scaling wrt the fiducial power spectrum (for var_NL_rescale=1.0 the code uses the Takahashi et al. (2012) version of halofit (see also Smith et al. 2003)).
 * 
 */

double ProjectedGalaxySample::compute_skewness_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale){
  return this->pointer_to_universe()->return_LOS_integrated_skewness(theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale, this->w_values, this->n_of_w_values);
}


/*
 * ProjectedGalaxySample::return_N_max_in_angular_tophat
 * 
 * Given a spherical volume of radius R at this->redshift, return the maximum N up to which the CiC histogram is computed in ProjectedGalaxySample::return_CiC_PDF. The python interface needs to know this before calling the CiC computation (which is why this function exists in the first place).
 * 
 */

int ProjectedGalaxySample::return_N_max_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale){
  
  double A = 2.0*constants::pi*(1.0-cos(theta_in_arcmin*constants::arcmin));
  double N_bar = A*this->density;
  double variance = this->compute_variance_in_angular_tophat(theta_in_arcmin, var_NL_rescale);
  
  return this->return_N_max(N_bar, variance);
  
}


/*
 * ProjectedGalaxySample::set_b2_to_minimise_negative_densities_in_angular_tophat
 * 
 * Changes quadratic_bias, trying to avoid negative values for the galaxy density. However, it also enforces the constraint that delta_g should be a monotonic function of delta_m.
 * 
 * 
 */

double ProjectedGalaxySample::set_b2_to_minimise_negative_densities_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale){
  
  double variance = this->compute_variance_in_angular_tophat(theta_in_arcmin, var_NL_rescale);
  return this->set_b2_to_minimise_negative_densities(variance);
  
}

/*
 * ProjectedGalaxySample::return_CiC_PDF_in_angular_tophat
 * 
 * Returns an array containing the probabilities of finding N galaxies in an angular tophat filter of radius theta.
 * 
 */

vector<double> ProjectedGalaxySample::return_CiC_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale){
  
  vector<vector<double> > PDF_data = this->pointer_to_universe()->compute_LOS_projected_PDF(this->w_values, this->n_of_w_values, theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale);
  
  double A = 2.0*constants::pi*(1.0-cos(theta_in_arcmin*constants::arcmin));
  double N_bar = A*this->density;
  
  return this->return_CIC_from_matter_density_PDF(N_bar, PDF_data);
  
}

/*
 * ProjectedGalaxySample::return_CiC_saddle_point_PDF_in_angular_tophat
 * 
 * Same as return_CiC_PDF_in_angular_tophat, but: to compute the PDF of matter density fluctuations, the CGF-to-PDF inverse Laplace transformation is approximated by it's saddle point. 
 * 
 */

vector<double> ProjectedGalaxySample::return_CiC_saddle_point_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale){
  
  vector<vector<double> > PDF_data = this->pointer_to_universe()->compute_LOS_projected_PDF_saddle_point(this->w_values, this->n_of_w_values, theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale);
  
  double A = 2.0*constants::pi*(1.0-cos(theta_in_arcmin*constants::arcmin));
  double N_bar = A*this->density;
  
  return this->return_CIC_from_matter_density_PDF(N_bar, PDF_data);
  
}


/*
 * ProjectedGalaxySample::return_joint_saddle_point_PDF_Ng_kappaCMB_in_angular_tophat
 * 
 * Return joint PDF of galaxy number counts and CMB convergence, smoothed over angular top-hat of radius theta_in_arcmin.
 * NOTE: this used a saddle point approximation for the joint PDF p(delta, kappa) instead of the full inverse Laplace transform.
 * 
 */

void ProjectedGalaxySample::return_joint_saddle_point_PDF_Ng_kappaCMB_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, vector<vector<double> > *Ng_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid){
  
  
  vector<vector<double> > d_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  this->pointer_to_universe()->compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point(theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale, kappa_min, kappa_max, this->w_values, this->n_of_w_values, &d_grid, &k_grid, &p_grid);

  int N_delta = d_grid.size(); 
  int N_kappa = d_grid[0].size(); // in principle N_kappa should be == N_delta. But let's be ignorant on what the other class is doing.
  
  double A = 2.0*constants::pi*(1.0-cos(theta_in_arcmin*constants::arcmin));
  double N_bar = A*this->density;
  double variance = this->compute_variance_in_angular_tophat(theta_in_arcmin, var_NL_rescale);
  double norm;
  double d_delta, integrand;
  int N_max = this->return_N_max(N_bar, variance);
  
  vector<vector<double> > P_of_N_given_delta(N_max+1, vector<double>(N_delta, 0.0));
  
  // ISSUE: you use here the knowledge, that all columns of d_grid are identical.
  for(int d = 0; d < N_delta; d++){
    norm = 0.0;
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] = this->return_P_of_N_given_delta(n, N_bar, d_grid[d][0], variance);
      norm += P_of_N_given_delta[n][d];
    }
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] /= norm;
    }
  }
  
  (*Ng_grid) = vector<vector<double> >(N_max+1, vector<double>(N_kappa, 0.0));
  (*kappa_grid) = vector<vector<double> >(N_max+1, vector<double>(N_kappa, 0.0));
  (*PDF_grid) = vector<vector<double> >(N_max+1, vector<double>(N_kappa, 0.0));
  
  for(int k = 0; k < N_kappa; k++){
    for(int n = 0; n < N_max+1; n++){
      (*Ng_grid)[n][k] = double(n);
      (*kappa_grid)[n][k] = k_grid[0][k];
      for(int d = 0; d < N_delta-1; d++){
        d_delta = d_grid[d+1][0]-d_grid[d][0];
        integrand = 0.5*(P_of_N_given_delta[n][d]*p_grid[d][k] + P_of_N_given_delta[n][d+1]*p_grid[d+1][k]);
        (*PDF_grid)[n][k] += d_delta*integrand;
      }
    }
  }
    
}

void ProjectedGalaxySample::return_joint_saddle_point_PDF_Ng_kappaCMB_noisy_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_CMB_noise_variance, vector<vector<double> > *Ng_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid){
  
  
  vector<vector<double> > d_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  this->pointer_to_universe()->compute_LOS_projected_PDF_incl_CMB_kappa_saddle_point(theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale, kappa_min, kappa_max, this->w_values, this->n_of_w_values, &d_grid, &k_grid, &p_grid);

  int N_delta = d_grid.size(); 
  int N_kappa = d_grid[0].size(); // in principle N_kappa should be == N_delta. But let's be ignorant on what the other class is doing.
  
  double A = 2.0*constants::pi*(1.0-cos(theta_in_arcmin*constants::arcmin));
  double N_bar = A*this->density;
  double variance = this->compute_variance_in_angular_tophat(theta_in_arcmin, var_NL_rescale);
  double norm;
  double d_delta, d_kappa, integrand;
  int N_max = this->return_N_max(N_bar, variance);
  
  vector<vector<double> > P_of_N_given_delta(N_max+1, vector<double>(N_delta, 0.0));
  
  // ISSUE: you use here the knowledge, that all columns of d_grid are identical.
  for(int d = 0; d < N_delta; d++){
    norm = 0.0;
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] = this->return_P_of_N_given_delta(n, N_bar, d_grid[d][0], variance);
      norm += P_of_N_given_delta[n][d];
    }
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] /= norm;
    }
  }
  
  (*Ng_grid) = vector<vector<double> >(N_max+1, vector<double>(N_kappa, 0.0));
  (*kappa_grid) = vector<vector<double> >(N_max+1, vector<double>(N_kappa, 0.0));
  vector<vector<double> > PDF_grid_no_noise(N_max+1, vector<double>(N_kappa, 0.0));
  (*PDF_grid) = vector<vector<double> >(N_max+1, vector<double>(N_kappa, 0.0));
  
  for(int k = 0; k < N_kappa; k++){
    for(int n = 0; n < N_max+1; n++){
      (*Ng_grid)[n][k] = double(n);
      (*kappa_grid)[n][k] = k_grid[0][k];
      for(int d = 0; d < N_delta-1; d++){
        d_delta = d_grid[d+1][0]-d_grid[d][0];
        integrand = 0.5*(P_of_N_given_delta[n][d]*p_grid[d][k] + P_of_N_given_delta[n][d+1]*p_grid[d+1][k]);
        PDF_grid_no_noise[n][k] += d_delta*integrand;
      }
    }
  }
  
  double kappa, kappa_tilde;
  double integrand_normalisation = 1.0/sqrt(2.0*constants::pi*kappa_CMB_noise_variance);
  for(int k = 0; k < N_kappa; k++){
    kappa = k_grid[0][k];
    for(int n = 0; n < N_max+1; n++){
      for(int kk = 0; kk < N_kappa-1; kk++){
        d_kappa = k_grid[0][kk+1]-k_grid[0][kk];
        kappa_tilde = k_grid[0][kk];
        integrand = 0.5*PDF_grid_no_noise[n][kk]*exp(-0.5*pow(kappa-kappa_tilde, 2)/kappa_CMB_noise_variance);
        kappa_tilde = k_grid[0][kk+1];
        integrand += 0.5*PDF_grid_no_noise[n][kk+1]*exp(-0.5*pow(kappa-kappa_tilde, 2)/kappa_CMB_noise_variance);
        integrand *= integrand_normalisation;
        (*PDF_grid)[n][k] += d_kappa*integrand;
      }
    }
  }
    
}

void ProjectedGalaxySample::return_joint_saddle_point_PDF_Ng_kappa_noisy_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale, double kappa_min, double kappa_max, double kappa_noise_variance, vector<double> w_values_lensing_kernel, vector<double> lensing_kernel_values, vector<vector<double> > *Ng_grid, vector<vector<double> > *kappa_grid, vector<vector<double> > *PDF_grid){
  
  
  vector<vector<double> > d_grid;
  vector<vector<double> > k_grid;
  vector<vector<double> > p_grid;
  
  vector<vector<double> > rebinned_data = return_joint_binning(this->w_values, this->n_of_w_values, w_values_lensing_kernel, lensing_kernel_values);
  
  for(int t = 0; t < rebinned_data[0].size(); t++){
    cout << t << "   ";
    cout << rebinned_data[0][t] << "   ";
    cout << rebinned_data[1][t] << "   ";
    cout << rebinned_data[2][t] << "\n";
  }
  
  this->pointer_to_universe()->compute_LOS_projected_PDF_incl_kappa_saddle_point(theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale, kappa_min, kappa_max, rebinned_data[0], rebinned_data[1], rebinned_data[2], &d_grid, &k_grid, &p_grid);

  int N_delta = d_grid.size(); 
  int N_kappa = d_grid[0].size(); // in principle N_kappa should be == N_delta. But let's be ignorant on what the other class is doing.
  
  double A = 2.0*constants::pi*(1.0-cos(theta_in_arcmin*constants::arcmin));
  double N_bar = A*this->density;
  double variance = this->compute_variance_in_angular_tophat(theta_in_arcmin, var_NL_rescale);
  double norm;
  double d_delta, d_kappa, integrand;
  int N_max = this->return_N_max(N_bar, variance);
  
  vector<vector<double> > P_of_N_given_delta(N_max+1, vector<double>(N_delta, 0.0));
  
  // ISSUE: you use here the knowledge, that all columns of d_grid are identical.
  for(int d = 0; d < N_delta; d++){
    norm = 0.0;
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] = this->return_P_of_N_given_delta(n, N_bar, d_grid[d][0], variance);
      norm += P_of_N_given_delta[n][d];
    }
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] /= norm;
    }
  }
  
  (*Ng_grid) = vector<vector<double> >(N_max+1, vector<double>(N_kappa, 0.0));
  (*kappa_grid) = vector<vector<double> >(N_max+1, vector<double>(N_kappa, 0.0));
  vector<vector<double> > PDF_grid_no_noise(N_max+1, vector<double>(N_kappa, 0.0));
  (*PDF_grid) = vector<vector<double> >(N_max+1, vector<double>(N_kappa, 0.0));
  
  for(int k = 0; k < N_kappa; k++){
    for(int n = 0; n < N_max+1; n++){
      (*Ng_grid)[n][k] = double(n);
      (*kappa_grid)[n][k] = k_grid[0][k];
      for(int d = 0; d < N_delta-1; d++){
        d_delta = d_grid[d+1][0]-d_grid[d][0];
        integrand = 0.5*(P_of_N_given_delta[n][d]*p_grid[d][k] + P_of_N_given_delta[n][d+1]*p_grid[d+1][k]);
        PDF_grid_no_noise[n][k] += d_delta*integrand;
      }
    }
  }
  
  if(kappa_noise_variance > 0.0){
    double kappa, kappa_tilde;
    double integrand_normalisation = 1.0/sqrt(2.0*constants::pi*kappa_noise_variance);
    for(int k = 0; k < N_kappa; k++){
      kappa = k_grid[0][k];
      for(int n = 0; n < N_max+1; n++){
        for(int kk = 0; kk < N_kappa-1; kk++){
          d_kappa = k_grid[0][kk+1]-k_grid[0][kk];
          kappa_tilde = k_grid[0][kk];
          integrand = 0.5*PDF_grid_no_noise[n][kk]*exp(-0.5*pow(kappa-kappa_tilde, 2)/kappa_noise_variance);
          kappa_tilde = k_grid[0][kk+1];
          integrand += 0.5*PDF_grid_no_noise[n][kk+1]*exp(-0.5*pow(kappa-kappa_tilde, 2)/kappa_noise_variance);
          integrand *= integrand_normalisation;
          (*PDF_grid)[n][k] += d_kappa*integrand;
        }
      }
    }
  }
  else{
    (*PDF_grid) = PDF_grid_no_noise;
  }
    
}




/*
 * ProjectedGalaxySample::return_kappa_PDF_in_angular_tophat
 * 
 * Returns an array whose 1st column contains values of convergence kappa and whose
 * 2nd column contains the PDF p(kappa).
 * 
 */

vector<vector<double> > ProjectedGalaxySample::return_kappa_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale){
  
  vector<vector<double> > PDF_data = this->pointer_to_universe()->compute_LOS_projected_PDF(this->w_values, this->lensing_kernel_values, theta_in_arcmin*constants::arcmin, f_NL, var_NL_rescale);
  
  return PDF_data;
  
}


/*
 * ProjectedGalaxySample::return_LOS_data
 * 
 * Returns bins in co-moving distance and the corresponding line-of-sight distribution of galaxies and lensing kernel.
 * 
 */
void ProjectedGalaxySample::return_LOS_data(vector<double> *w_vals, vector<double> *n_of_w_vals, vector<double> *lensing_kernel_vals){
  (*w_vals) = this->w_values;
  (*n_of_w_vals) = this->n_of_w_values;
  (*lensing_kernel_vals) = this->lensing_kernel_values;
}


