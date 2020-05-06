#include "GalaxySample.h"


GalaxySample::GalaxySample(Matter* matter, double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1){
  
  this->matter = matter;
  this->universe = this->matter->universe;
  this->cosmology = this->universe->return_cosmology(); 
  
  this->density = density_in_Mpc_over_h_cubed*pow(c_over_e5,3); // changing to units c/H_0 == 1
  this->redshift = z;
  this->linear_bias = b1;
  this->quadratic_bias = b2;
  this->alpha_0 = a0;
  this->alpha_1 = a1;
  
}


GalaxySample::~GalaxySample(){
}


/*
 * GalaxySample::change_parameters
 * 
 * Changes the parameters that characterise the galaxy sample.
 * 
 * 
 */

void GalaxySample::change_parameters(double z, double density_in_Mpc_over_h_cubed, double b1, double b2, double a0, double a1){
  this->density = density_in_Mpc_over_h_cubed*pow(c_over_e5,3); // changing to units c/H_0 == 1
  this->redshift = z;
  this->linear_bias = b1;
  this->quadratic_bias = b2;
  this->alpha_0 = a0;
  this->alpha_1 = a1;
}


/*
 * GalaxySample::set_b2_to_minimise_negative_densities
 * 
 * Changes quadratic_bias, trying to avoid negative values for the galaxy density. However, it also enforces the constraint that delta_g should be a monotonic function of delta_m.
 * 
 * 
 */

double GalaxySample::set_b2_to_minimise_negative_densities(double z, double R_in_Mpc_over_h, double var_NL_rescale){
  
  double var = var_NL_rescale*this->matter->return_non_linear_variance(z, R_in_Mpc_over_h);
  double R = R_in_Mpc_over_h/c_over_e5;
  
  if(this->linear_bias>1.0)
    this->quadratic_bias = min(0.5*this->linear_bias,(this->linear_bias-1.0)/(1.0-var));
  else
    this->quadratic_bias = 0.0;
  
  cout << " b1 = " << this->linear_bias << '\n';
  cout << " b2 = " << this->quadratic_bias << '\n';
  cout << "var = " << var << '\n';
  
  return this->quadratic_bias;
}


/*
 * GalaxySample::P_of_N_given_delta
 * 
 * Given the matter density contrast delta within a volume V, compute the probability of finding N galaxies in that volume. The variance is only relevant when quadratic_bias != 0.
 * 
 */

double GalaxySample::return_P_of_N_given_delta(int N, double V, double delta, double variance){
  double delta_g = this->linear_bias*delta + this->quadratic_bias*(delta*delta - variance);
  if(delta_g<-1){
    delta_g=-1.0;
    if(this->error_flag_negative_density==0){
      this->error_flag_negative_density = 1;
      cerr << "delta_g found to be < -1 in GalaxySample::P_of_N_given_delta!\n";
      cerr << "Set to delta_g == -1.\n";
    }
  }
  
  if(delta_g<=-1.) {
    if(N > 0) return 0.;
    return 1.;
  }
  
  double galaxy_per_Poisson_halo = this->alpha_0 + this->alpha_1*delta; // In the convention of https://arxiv.org/abs/1710.05162 this uses indeed delta==delta_matter.
  double N_bar = this->density*V;
  
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
 * GalaxySample::return_N_max
 * 
 * Given a spherical volume of radius R at redshift z, return the maximum N up to which the CiC histogram is computed in GalaxySample::return_CiC_PDF. The python interface needs to know this before calling the CiC computation (which is why this function exists in the first place).
 * 
 */

int GalaxySample::return_N_max(double z, double R_in_Mpc_over_h, double var_NL_rescale){
  
  double var = var_NL_rescale*this->matter->return_non_linear_variance(z, R_in_Mpc_over_h);
  double R = R_in_Mpc_over_h/c_over_e5;
  
  // to estimate the maximum N:
  // assume delta_matter is log-normal with delta_0 = -1
  // go to 5\sigma of both delta and shot-noise
  
  double var_Gauss = log(1.0+var);
  double delta_max = (exp(-0.5*var_Gauss + 5.0*sqrt(var_Gauss))-1.0);
  double delta_g_max = this->linear_bias*delta_max + this->quadratic_bias*(delta_max*delta_max - var);
  double N_bar = 4.0*constants::pi/3.0*pow(R,3)*delta_g_max*this->density;
  double galaxies_per_Poisson_halo = this->alpha_0 + delta_max*this->alpha_1;
  
  return int(0.5+N_bar+5.0*sqrt(N_bar*galaxies_per_Poisson_halo));
  
}

int GalaxySample::return_N_max_and_variance(double z, double R_in_Mpc_over_h, double* variance){
  
  (*variance) = this->matter->return_non_linear_variance(z, R_in_Mpc_over_h);
  double R = R_in_Mpc_over_h/c_over_e5;
  
  // to estimate the maximum N:
  // assume delta_matter is log-normal with delta_0 = -1
  // go to 5\sigma of both delta and shot-noise
  
  double var_Gauss = log(1.0+(*variance));
  double delta_max = (exp(-0.5*var_Gauss + 5.0*sqrt(var_Gauss))-1.0);
  double delta_g_max = this->linear_bias*delta_max + this->quadratic_bias*(delta_max*delta_max - (*variance));
  double N_bar = 4.0*constants::pi/3.0*pow(R,3)*delta_g_max*this->density;
  double galaxies_per_Poisson_halo = this->alpha_0 + delta_max*this->alpha_1;
  
  return int(0.5+N_bar+5.0*sqrt(N_bar*galaxies_per_Poisson_halo));
  
}





/*
 * GalaxySample::return_CiC_PDF
 * 
 * Returns an array containing the probabilities of finding N galaxies in a spherical colume of radius R_in_Mpc_over_h at redshift z.
 * 
 */

vector<double> GalaxySample::return_CiC_PDF(double z, double R_in_Mpc_over_h, double f_NL, double var_NL_rescale){
  
  vector<vector<double> > PDF_data = this->matter->compute_PDF_3D(z, R_in_Mpc_over_h, f_NL, var_NL_rescale);
  
  double V = 4.0*constants::pi/3.0*pow(R_in_Mpc_over_h/c_over_e5,3);
  double d_delta = PDF_data[0][1] - PDF_data[0][0];
  double variance;
  
  int N_max = this->return_N_max_and_variance(z, R_in_Mpc_over_h, &variance);
  variance *= var_NL_rescale;
  int n_delta = PDF_data[0].size();
  
  vector<double> P_of_N(N_max+1, 0.0);
  vector<vector<double> > P_of_N_given_delta(N_max+1, vector<double>(n_delta, 0.0));
  
  for(int d = 0; d < n_delta; d++){
    double norm = 0.0;
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] = this->return_P_of_N_given_delta(n, V, PDF_data[0][d], variance);
      norm += P_of_N_given_delta[n][d];
    }
    for(int n = 0; n < N_max+1; n++){
      P_of_N_given_delta[n][d] /= norm;
    }
  }
  
  for(int n = 0; n < N_max+1; n++){
    for(int d = 0; d < n_delta; d++){
      P_of_N[n] += P_of_N_given_delta[n][d]*PDF_data[1][d];
    }
    P_of_N[n] *= d_delta;
  }
  
  return P_of_N;
  
}












