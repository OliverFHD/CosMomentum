

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
  
  double a_1, a_2;
  double w_1, w_2;
  double w_final = this->w_values[Nz-1];
  
  // now calculating the lensing kernel (if ProjectedGalaxySample is used as source galaxy sample)
  for(int i = 0; i < Nz-1; i++){
    w = 0.5*(this->w_values[i]+this->w_values[i+1]);
    w_1 = this->w_values[i+1];
    dw = this->w_values[i+1] - w;
    a_1 = this->pointer_to_universe()->a_at_eta(eta_0 - w_1);
    this->lensing_kernel_values[i] += dw*0.5*(w*(w_1-w)/a_1/w_1);
    for(int j = i+1; j < Nz-1; j++){
      w_1 = this->w_values[j];
      w_2 = this->w_values[j+1];
      dw = w_2 - w_1;
      a_1 = this->pointer_to_universe()->a_at_eta(eta_0 - w_1);
      a_2 = this->pointer_to_universe()->a_at_eta(eta_0 - w_2);
      //w*(w_final-w)/w_final/a
      this->lensing_kernel_values[i] += dw*0.5*(w*(w_1-w)/a_1/w_1);
      this->lensing_kernel_values[i] += dw*0.5*(w*(w_2-w)/a_2/w_2);
    }
    this->lensing_kernel_values[i] *= 1.5*this->pointer_to_universe()->return_Omega_m();
  
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





