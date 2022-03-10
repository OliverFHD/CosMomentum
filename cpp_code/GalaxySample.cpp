
GalaxySample::GalaxySample(FlatInhomogeneousUniverseLCDM* universe, double b1, double b2, double a0, double a1, BIAS_MODEL b_model){

  this->set_universe(universe);
  this->set_parameters(b1, b2, a0, a1, b_model);

  if(this->bias_model == LAGRANGIAN){
    cout << "LAGRANGIAN bias model.\n";
  }
  else{
    cout << "EULERIAN bias model.\n";
  }

  N_samples=0;
  delta_m_samples=0;
  b_samples=1;
  n_samples=0;

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

void GalaxySample::set_parameters(double b1, double b2, double a0, double a1, BIAS_MODEL b_model){
  this->bias_model = b_model;
  if(b_model == EULERIAN){
    this->linear_bias = b1;
    this->quadratic_bias = b2;
    this->linear_Lagrangian_bias = -1.0;
    this->quadratic_Lagrangian_bias = -1.0;
  }
  else{
    this->linear_Lagrangian_bias = b1;
    this->quadratic_Lagrangian_bias = b2;
    this->linear_bias = -1.0;
    this->quadratic_bias = -1.0;
  }
  this->alpha_0 = a0;
  this->alpha_1 = a1;
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

  double b_lin = this->linear_bias;
  double b_quad = this->quadratic_bias;
  if(this->bias_model == LAGRANGIAN){
    b_lin = 1.0 + this->linear_Lagrangian_bias;
    b_quad = this->quadratic_Lagrangian_bias + 8.0/21.0*this->linear_Lagrangian_bias;
  }

  double b_quad_new = 2.0*min(0.5*b_lin,(b_lin-1.0)/(1.0-variance));

  if(this->bias_model == LAGRANGIAN){
    this->quadratic_Lagrangian_bias += b_quad_new - b_quad;
    return this->quadratic_Lagrangian_bias;
  }
  else{
    this->quadratic_bias = b_quad_new;
  }

  return b_quad_new;
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
  double b_lin = this->linear_bias;
  double b_quad = this->quadratic_bias;
  if(this->bias_model == LAGRANGIAN){
    // Approximate relations between b_Euler and b_Lagrange
    // (Not vital that this is 100% accurate, since we're only
    // calculating N_max up to which CiC histogram is calculated.)
    b_lin = 1.0 + this->linear_Lagrangian_bias;
    b_quad = this->quadratic_Lagrangian_bias + 8.0/21.0*this->linear_Lagrangian_bias;
  }
  double delta_g_max = b_lin*delta_max + b_quad/2.0*(delta_max*delta_max - variance);
  if(b_quad < 0.0){
    delta_g_max = -0.5*b_lin/(b_quad/2.0);
  }

  double N_bar_max = N_bar*(1.0+delta_g_max);
  double galaxies_per_Poisson_halo = this->alpha_0 + delta_max*this->alpha_1;

  return int(0.5+N_bar_max+5.0*sqrt(N_bar_max*galaxies_per_Poisson_halo));

}


/*
 * GalaxySample::return_bias_model
 *
 * Returns GalaxySample.return_bias_model (= EULERIAN or LAGRANGIAN).
 *
 */

BIAS_MODEL GalaxySample::return_bias_model(){
  return this->bias_model;
};



/*
 * GalaxySample::delta_g_Eulerian
 *
 * Given the matter density contrast delta within a volume V, return the galaxy density contrast in the same volume.
 * The variance is meant to be <\delta_{m,V}^2> and is only relevant when quadratic_bias != 0.
 * This function uses the Eulerian bias parametrisation.
 *
 */

double GalaxySample::delta_g_Eulerian(double delta, double variance){
  return this->linear_bias*delta + this->quadratic_bias/2.0*(delta*delta - variance);
}


/*
 * GalaxySample::delta_g_Lagrangian
 *
 * Given the matter density contrast delta within a volume V, return the galaxy density contrast in the same volume.
 * PDF_data is in the PDF-output format of class FlatInhomogeneousUniverseLCDM.
 * This function uses the Lagrangian bias parametrisation.
 *
 */

vector<vector<double> > GalaxySample::delta_g_Lagrangian(vector<vector<double> > *PDF_data){
  int N_delta = (*PDF_data)[0].size();
  vector<vector<double> > delta_m_and_g(2, vector<double>(N_delta, 0.0));
  delta_m_and_g[0] = (*PDF_data)[0];

  for(int d = 0; d < N_delta; d++){
    delta_m_and_g[1][d] = (*PDF_data)[0][d] + (*PDF_data)[2][d]*this->linear_Lagrangian_bias + (*PDF_data)[3][d]*this->quadratic_Lagrangian_bias/2.0;
  }
  return delta_m_and_g;
}


/*
 * GalaxySample::return_P_of_N_given_delta_g
 *
 * Given the (shot-noise free) galaxy density contrast delta_g within a volume V, compute the probability of finding N galaxies in that volume.
 *
 */

double GalaxySample::return_P_of_N_given_delta_g(int N, double N_bar, double delta_m, double delta_g){

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

  if(N_samples) { // samples of N,delta have been stored - use them!
    return_P_of_N_given_delta_g_sampled(N, delta_g);
  }

  double galaxy_per_Poisson_halo = this->alpha_0 + this->alpha_1*delta_m; // In the convention of https://arxiv.org/abs/1710.05162 this uses indeed delta==delta_matter.

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

/*
 * GalaxySample::store_N_delta_samples
 *
 * Store a given list of samples of N and delta_m in e.g. a mock, and an associated linear bias,
 * and in any future call of return_P_of_N_given_delta_g actually call return_P_of_N_given_delta_g_sampled on these
 *
 */
void GalaxySample::store_N_delta_samples(double db_samples, int dn_samples, double *dN_samples, double *ddelta_m_samples) {

  b_samples = db_samples;
  n_samples = dn_samples;
  N_samples = dN_samples;
  delta_m_samples = ddelta_m_samples;

}

/*
 * GalaxySample::return_P_of_N_given_delta_g_sampled
 *
 * Given a list of samples of N and delta_m in e.g. a mock, and an associated linear bias,
 * return the probability of finding N galaxies in a volume with (shot-noise free) galaxy density contrast delta_g.
 * Caveat: does not ensure the mean of N is what you would like it to be (but does ensure linear bias is).
 *
 */

double GalaxySample::return_P_of_N_given_delta_g_sampled(int N, double delta_g){

  double delta_m = delta_g/b_samples; // this is the delta_m at which we would like to evaluate the sampled PDF

  // determine bandwidth

  double width_N = pow(n_samples,-1./6)*gsl_stats_sd(N_samples, 1, n_samples); // Silverman's rule of thumb in 2D
  double width_delta = pow(n_samples,-1./6)*gsl_stats_sd(delta_m_samples, 1, n_samples);

  // Epanechnikov kernel
  double p = 0.;
  for(int i=0; i<n_samples; i++) {
    double dN = fabs(N-N_samples[i]);
    if(dN>=width_N) continue;
    double ddelta = fabs(delta_m - delta_m_samples[i]);
    double x = pow(dN/width_N,2)+pow(ddelta/width_delta,2);
    if(x>1) continue;
    p += (1.-x*x);
  }

  return p;
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

  /*
   * If using Lagrangian bias parametrisation,
   * then we need information on bias_term_1 and bias_term_2
   * (cf. class FlatInhomogeneousUniverseLCDM)
   *
   */
  vector<vector<double> > delta_m_and_g;
  if(this->return_bias_model() == LAGRANGIAN){
    delta_m_and_g = this->delta_g_Lagrangian(&PDF_data);
  }

  vector<double> P_of_N(N_max+1, 0.0);
  vector<vector<double> > P_of_N_given_delta(N_max+1, vector<double>(n_delta, 0.0));

  double norm;
  double delta_g;
  for(int d = 0; d < n_delta; d++){
    norm = 0.0;
    if(this->return_bias_model() == EULERIAN){
      delta_g = this->delta_g_Eulerian(PDF_data[0][d], variance);
    }
    else if(this->return_bias_model() == LAGRANGIAN){
      delta_g = interpolate_neville_aitken(PDF_data[0][d], &delta_m_and_g[0], &delta_m_and_g[1], constants::order_of_interpolation);
    }
    else{
      error_handling::general_error_message("Invalid value for GalaxySample.bias_model in function GalaxySample.return_CIC_from_matter_density_PDF .");
    }

    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] = this->return_P_of_N_given_delta_g(n, N_bar, PDF_data[0][d], delta_g);
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
