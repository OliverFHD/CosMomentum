

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
 * Read in the redshift distribution of the galaxy sample and turn it into a co-moving distance distribution. The textfile "n_of_z_input_file" should have two columns, where column 1 contains the lower-redshift edge of each histogram bin and column 2 contains the redshift histrogram (the histogram doesn't need to be normalised, since normalisation is enforced later on in the code).
 * 
 */
void ProjectedGalaxySample::set_n_of_w_data(string n_of_z_input_file){
  
  vector<vector<double> > nofz_data = read_table_double(n_of_z_input_file);
  
  int Nz = nofz_data.size();
  this->w_values = vector<double>(Nz, 0.0);
  this->n_of_w_values = vector<double>(Nz, 0.0);
  double eta_0 = this->pointer_to_universe()->eta_at_a(1.0);
  double scale, dz, dw;
  
  for(int i = 0; i < Nz; i++){
    scale = 1.0/(1.0+nofz_data[i][0]);
    this->w_values[i] = eta_0 - this->pointer_to_universe()->eta_at_a(scale);
  }
  
  double norm = 0.0;
  for(int i = 0; i < Nz-1; i++){
    dz = nofz_data[i+1][0]-nofz_data[i][0];
    dw = this->w_values[i+1]-this->w_values[i];
    this->n_of_w_values[i] = nofz_data[i][1]*dz/dw;
    norm += nofz_data[i][1]*dz;
    cout << i <<  "  ";
    cout << this->w_values[i] <<  "  ";
    cout << this->n_of_w_values[i] <<  "\n";
  }
  
  for(int i = 0; i < Nz-1; i++) this->n_of_w_values[i] /= norm;
  
}


/*
 * ProjectedGalaxySample::compute_variance_in_angular_tophat
 * 
 * Compute variance of the projected density contrast smoothed with a angular top-hat filter of radius theta. The variance allows for a re-scaling wrt the fiducial power spectrum (for var_NL_rescale=1.0 the code uses the Takahashi et al. (2012) version of halofit (see also Smith et al. 2003)).
 * 
 */

double ProjectedGalaxySample::compute_variance_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale){
  return var_NL_rescale*this->pointer_to_universe()->return_LOS_integrated_variance(theta_in_arcmin*constants::arcmin, this->w_values, this->n_of_w_values);
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



