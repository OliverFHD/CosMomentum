

GalaxySample3D::GalaxySample3D(FlatInhomogeneousUniverseLCDM* universe, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1, BIAS_MODEL b_model) : GalaxySample(universe, b1, b2, a0, a1, b_model){
  
  this->redshift = z;
  this->density = density_in_Mpc_over_h_cubed*pow(constants::c_over_e5, 3); // changing from Mpc/h to units with c/H_0 == 1
  
}


GalaxySample3D::~GalaxySample3D(){
}

void GalaxySample3D::set_parameters_3D(double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1, BIAS_MODEL b_model){
  this->density = density_in_Mpc_over_h_cubed*pow(constants::c_over_e5, 3); // changing from Mpc/h to units with c/H_0 == 1
  this->redshift = z;
  this->set_parameters(b1, b2, a0, a1, b_model);
}


/*
 * GalaxySample3D::compute_variance_in_3D_tophat
 * 
 * Compute variance of the density contrast at this->redshift smoothed with a spherical top-hat filter of radius R_in_Mpc_over_h. The variance allows for a re-scaling wrt the fiducial power spectrum (for var_NL_rescale=1.0 the code uses the Takahashi et al. (2012) version of halofit (see also Smith et al. 2003)).
 * 
 */

double GalaxySample3D::compute_variance_in_3D_tophat(double R_in_Mpc_over_h, double var_NL_rescale){
  return this->pointer_to_universe()->return_non_linear_variance(this->redshift, R_in_Mpc_over_h, var_NL_rescale);
}


/*
 * GalaxySample3D::return_N_max_in_3D_tophat
 * 
 * Given a spherical volume of radius R at this->redshift, return the maximum N up to which the CiC histogram is computed in GalaxySample3D::return_CiC_PDF. The python interface needs to know this before calling the CiC computation (which is why this function exists in the first place).
 * 
 */

int GalaxySample3D::return_N_max_in_3D_tophat(double R_in_Mpc_over_h, double var_NL_rescale){
  
  double V = 4.0*constants::pi/3.0*pow(R_in_Mpc_over_h/constants::c_over_e5, 3);
  double N_bar = V*this->density;
  double variance = this->compute_variance_in_3D_tophat(R_in_Mpc_over_h, var_NL_rescale);
  
  return this->return_N_max(N_bar, variance);
  
}


/*
 * GalaxySample3D::set_b2_to_minimise_negative_densities_in_3D_tophat
 * 
 * Changes quadratic_bias, trying to avoid negative values for the galaxy density. However, it also enforces the constraint that delta_g should be a monotonic function of delta_m.
 * 
 * 
 */

double GalaxySample3D::set_b2_to_minimise_negative_densities_in_3D_tophat(double R_in_Mpc_over_h, double var_NL_rescale){
  
  double variance = this->compute_variance_in_3D_tophat(R_in_Mpc_over_h, var_NL_rescale);
  return this->set_b2_to_minimise_negative_densities(variance);
  
}

/*
 * GalaxySample3D::return_CiC_PDF_in_3D_tophat
 * 
 * Returns an array containing the probabilities of finding N galaxies in a spherical colume of radius R_in_Mpc_over_h at this->redshift.
 * 
 */

vector<double> GalaxySample3D::return_CiC_PDF_in_3D_tophat(double R_in_Mpc_over_h, double f_NL, double var_NL_rescale){
  
  vector<vector<double> > PDF_data = this->pointer_to_universe()->compute_PDF_3D(this->redshift, R_in_Mpc_over_h, f_NL, var_NL_rescale);
  
  double V = 4.0*constants::pi/3.0*pow(R_in_Mpc_over_h/constants::c_over_e5, 3);
  double N_bar = V*this->density;
  
  return this->return_CIC_from_matter_density_PDF(N_bar, PDF_data);
  
}












