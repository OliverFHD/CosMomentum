
GalaxySample::GalaxySample(FlatInhomogeneousUniverseLCDM* universe, double b1, double b2, double a0, double a1){
  
  this->set_universe(universe);
  this->set_parameters(b1, b2, a0, a1);
  
}


GalaxySample::~GalaxySample(){
}


/*
 * GalaxySample::set_parameters
 * 
 * Set the parameters that characterise the galaxy sample.
 * 
 * 
 */

void GalaxySample::set_parameters(double b1, double b2, double a0, double a1){
  this->linear_bias = b1;
  this->quadratic_bias = b2;
  this->alpha_0 = a0;
  this->alpha_1 = a1;
}


/*
 * GalaxySample::set_bias_model_from_br_parametrisation
 * 
 * Translating the br-parametrisation of bias and shot-noise (cf. https://arxiv.org/pdf/1710.05162.pdf) to the b-alpha_0-alpha_1 parametrisation.
 * 
 * 
 */

void GalaxySample::set_bias_model_from_br_parametrisation(double b_tilde, double r, double N_bar, double variance, double skewness){
  
  this->linear_bias = b_tilde*r;
  //if(this->quadratic_bias != 0.0){
  //  this->quadratic_bias = 0.0;
  //  cerr << "CAREFUL: currently, b_2 = 0 is still enforced when calling set_bias_model_from_br_parametrisation.\n";
  //}
  double delta_m0 = lognormal_tools::get_delta0(variance, skewness);
  
  this->alpha_0 = lognormal_tools::return_alpha_0(r, b_tilde, N_bar, variance, delta_m0);
  this->alpha_1 = lognormal_tools::return_alpha_1(r, b_tilde, N_bar, variance, delta_m0);
  cout << this->alpha_0 << '\n';
  cout << this->alpha_1 << '\n';
  cout << variance << '\n';
  cout << skewness << '\n';
  cout << delta_m0 << '\n';
}


/*
 * GalaxySample::set_universe
 * 
 * Anchor the galaxy sample in a universe.
 * 
 * 
 */

void GalaxySample::set_universe(FlatInhomogeneousUniverseLCDM* universe){
  this->universe = universe;
  this->cosmology = this->universe->return_cosmology(); 
}


/*
 * GalaxySample::set_b2_to_minimise_negative_densities
 * 
 * Changes quadratic_bias, trying to avoid negative values for the galaxy density. However, it also enforces the constraint that delta_g should be a monotonic function of delta_m.
 * 
 * 
 */

double GalaxySample::set_b2_to_minimise_negative_densities(double variance){
  
  if(this->linear_bias>1.0)
    this->quadratic_bias = min(0.5*this->linear_bias,(this->linear_bias-1.0)/(1.0-variance));
  else
    this->quadratic_bias = 0.0;
  
  cout << " b1 = " << this->linear_bias << '\n';
  cout << " b2 = " << this->quadratic_bias << '\n';
  cout << "var = " << variance << '\n';
  
  return this->quadratic_bias;
}


/*
 * GalaxySample::return_N_max
 * 
 * Given a spherical volume of radius R at this->redshift, return the maximum N up to which the CiC histogram is computed in GalaxySample::return_CiC_PDF. The python interface needs to know this before calling the CiC computation (which is why this function exists in the first place).
 * 
 */

int GalaxySample::return_N_max(double N_bar, double variance){
  
  // to estimate the maximum N:
  // assume delta_matter is log-normal with delta_0 = -1
  // go to 5\sigma of both delta and shot-noise
  
  double var_Gauss = log(1.0+variance);
  double delta_max = (exp(-0.5*var_Gauss + 5.0*sqrt(var_Gauss))-1.0);
  double delta_g_max = this->linear_bias*delta_max + this->quadratic_bias*(delta_max*delta_max - variance);
  if(this->quadratic_bias < 0.0){
    delta_g_max = -0.5*this->linear_bias/this->quadratic_bias;
  }
  
  double N_bar_max = N_bar*(1.0+delta_g_max);
  double galaxies_per_Poisson_halo = this->alpha_0 + delta_max*this->alpha_1;
  
  return int(0.5+N_bar_max+5.0*sqrt(N_bar_max*galaxies_per_Poisson_halo));
  
}



/*
 * GalaxySample::P_of_N_given_delta
 * 
 * Given the matter density contrast delta within a volume V, compute the probability of finding N galaxies in that volume. The variance is only relevant when quadratic_bias != 0.
 * 
 */

double GalaxySample::return_P_of_N_given_delta(int N, double N_bar, double delta, double variance){
  double delta_g = this->linear_bias*delta + this->quadratic_bias*(delta*delta - variance);
  
  if(delta_g<-1){
    delta_g=-1.0;
    if(this->get_error_flag_negative_density()==0){
      this->set_error_flag_negative_density(1);
      cerr << "delta_g found to be < -1 in GalaxySample::P_of_N_given_delta!\n";
      cerr << "Set to delta_g == -1.\n";
    }
  }
  
  if(delta_g<=-1.) {
    if(N > 0) return 0.;
    return 1.;
  }
  
  double galaxy_per_Poisson_halo = this->alpha_0 + this->alpha_1*delta; // In the convention of https://arxiv.org/abs/1710.05162 this uses indeed delta==delta_matter.
  
  if(galaxy_per_Poisson_halo == 1.0) return gsl_ran_poisson_pdf(N, N_bar*(1.0+delta_g)); // galaxy_per_Poisson_halo == 1.0 ==> shot-noise is Poisson
  if(galaxy_per_Poisson_halo <= 0.0) {
    // (unphysical) limit corresponding to no shot-noise
    if(int(N_bar*(1.0+delta_g)+0.5)==N) return 1.;
    return 0.0;
  }
  
  double x = double(N)/galaxy_per_Poisson_halo; // treat x as a Poisson variable
  double lambda = N_bar/galaxy_per_Poisson_halo*(1.0+delta_g);
  
  // --> ISSUE: this PDF is only perfectly normalised for integer galaxy_per_Poisson_halo. As of now one needs to remember this and re-normalise the entire CiC PDF once it has been computed. (Deviations from normalisation are expected to be small though.)
  return exp(x*log(lambda) - lambda - gsl_sf_lngamma(x + 1.0))/galaxy_per_Poisson_halo;
}



vector<double> GalaxySample::return_CIC_from_matter_density_PDF(double N_bar, vector<vector<double> > PDF_data){
    
  int n_delta = PDF_data[0].size();
  double d_delta;
  double integrand;
  double variance = 0.0;
  
  for(int d = 0; d < n_delta-1; d++){
    d_delta = PDF_data[0][d+1]-PDF_data[0][d];
    integrand = 0.5*(pow(PDF_data[0][d],2)*PDF_data[1][d]+pow(PDF_data[0][d+1],2)*PDF_data[1][d+1]);
    variance += integrand*d_delta;
  }
  
  // ISSUE: maybe make variance a parameter of this function instead
  int N_max = this->return_N_max(N_bar, variance);
  
  
  vector<double> P_of_N(N_max+1, 0.0);
  vector<vector<double> > P_of_N_given_delta(N_max+1, vector<double>(n_delta, 0.0));
  
  for(int d = 0; d < n_delta; d++){
    double norm = 0.0;
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] = this->return_P_of_N_given_delta(n, N_bar, PDF_data[0][d], variance);
      norm += P_of_N_given_delta[n][d];
    }
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] /= norm;
    }
  }
  
  for(int n = 0; n < N_max+1; n++){
    for(int d = 0; d < n_delta-1; d++){
      d_delta = PDF_data[0][d+1]-PDF_data[0][d];
      integrand = 0.5*(P_of_N_given_delta[n][d]*PDF_data[1][d] + P_of_N_given_delta[n][d+1]*PDF_data[1][d+1]);
      P_of_N[n] += d_delta*integrand;
    }
  }
  
  return P_of_N;
  
}







