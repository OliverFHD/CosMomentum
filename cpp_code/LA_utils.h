#include <gsl/gsl_fft_complex.h>

using namespace std;




/*
 * 
 * 
 *   \int_{-xMin}^{xMax} dx/(2\pi) f(x) e^{-ikx}
 * 
 * now choose k = 2\pi l/(n*\DeltaX)
 * 
 * = \sum_{j=-n/2}^{n/2-1} \DeltaX/(2\pi) f_j e^{-2\pi i j l/n}
 * = \sum_{j=0}^{n-1} \DeltaX/(2\pi) f_{j-n/2} e^{-2\pi i (j-n/2) l/n}
 * = e^{\pi i l}\DeltaX/(2\pi) \sum_{j=0}^{n-1} f_{j-n/2} e^{-2\pi i j l/n}
 * 
 * to do your Fourier, use
 * e^{\pi i l}\DeltaX/(2\pi) \sum_{j=0}^{n-1} f_{j-n/2} e^{\pi i j } e^{-2\pi i j l/n}
 * 
 */



extern "C" void fourier_transform(double* x, double* f_re, double* f_im, double* k, double* f_tilde_re, double* f_tilde_im, int log2n){
  
  int n = pow(2, log2n);
  int n_half = pow(2, log2n-1);
  double data[2*n];
  double prefactor;
  complex<double> dummy;
  complex<double> imag(0.0, 1.0);
  
  for(int i = 0; i < n; i++){
    data[2*i] = f_re[i]*pow(-1, i);
    data[2*i+1] = f_im[i]*pow(-1, i);
  }
  
  gsl_fft_complex_radix2_forward(data, 1, n);
  
  double dx = x[1]-x[0];
  double k_width = 2.0*constants::pi/dx;
  double k_min = -k_width/2.0;
  double dk = k_width/n;
  double dx_over_2pi = dx/(2.0*constants::pi);
  
  for(int i = 0; i < n; i++){
    prefactor = pow(-1, i - n_half)*dx_over_2pi;
    f_tilde_re[i] = data[2*i]*prefactor;
    f_tilde_im[i] = data[2*i+1]*prefactor;
  }
  
  for(int i = 0; i < n; i++){
    k[i] = k_min + i*dk;
  }
  
}


extern "C" void fourier_transform_2D(double* x, double* y, double* F_re, double* F_im, double* kx, double* ky, double* F_tilde_re, double* F_tilde_im, int log2n){
  
  int n = pow(2, log2n);
  int n_half = pow(2, log2n-1);
  double dummy_re[n];
  double dummy_im[n];
  double dummy_tilde_re[n];
  double dummy_tilde_im[n];
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      dummy_re[j] = F_re[i*n+j];
      dummy_im[j] = F_im[i*n+j];
    }
    // fourier_transform(lambda_values, integrand_re, integrand_im, d_values,           p_re,           p_im, log2n);
    fourier_transform   (x            ,     dummy_re,     dummy_im,       kx, dummy_tilde_re, dummy_tilde_im, log2n);
    for(int j = 0; j < n; j++){
      F_re[i*n+j] = dummy_tilde_re[j];
      F_im[i*n+j] = dummy_tilde_im[j];
    }
  }
  
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      dummy_re[j] = F_re[j*n+i];
      dummy_im[j] = F_im[j*n+i];
    }
    fourier_transform(y, dummy_re, dummy_im, ky, dummy_tilde_re, dummy_tilde_im, log2n);
    for(int j = 0; j < n; j++){
      F_tilde_re[j*n+i] = dummy_tilde_re[j];
      F_tilde_im[j*n+i] = dummy_tilde_im[j];
    }
  }
  
}

/*
 * return_joint_binning
 * 
 * - Given 2 binned historgrams kernel_1(z_vals_1) and kernel_2(z_vals_2), where z_vals_* are
 *   the histogram bin edges, this function rebins both histograms so that they live on the
 *   same bins.
 * 
 */
vector<vector<double> > return_joint_binning(vector<double> z_vals_1, vector<double> kernel_1, vector<double> z_vals_2, vector<double> kernel_2){
  
  vector<double> z_vals_new(0, 0.0);
  double current_z = min(z_vals_1[0], z_vals_2[0]);
  int index_1 = 0;
  int index_2 = 0;
  int Nz_1 = z_vals_1.size();
  int Nz_2 = z_vals_2.size();
  
  while(index_1 < Nz_1 && index_2 < Nz_2){
    if(z_vals_1[index_1] == z_vals_2[index_2]){
      current_z = z_vals_1[index_1];
      index_1++;
      index_2++;
    }
    else if(z_vals_1[index_1] < z_vals_2[index_2]){
      current_z = z_vals_1[index_1];
      index_1++;
    }
    else{
      current_z = z_vals_2[index_2];
      index_2++;
    }
    z_vals_new.push_back(current_z);
  }
  
  while(index_1 < Nz_1){
    current_z = z_vals_1[index_1];
    z_vals_new.push_back(current_z);
    index_1++;
  }
  
  while(index_2 < Nz_2){
    current_z = z_vals_2[index_2];
    z_vals_new.push_back(current_z);
    index_2++;
  }
  
  int Nz_new = z_vals_new.size();
  // 1st row: new bin edges
  // 2nd row: rebinned kernel_1
  // 3rd row: rebinned kernel_2
  vector<vector<double> > rebinned_data(3, vector<double>(Nz_new, 0.0));
  rebinned_data[0] = z_vals_new;
  
  for(int z = 0; z < Nz_new; z++){
    current_z = z_vals_new[z];
    index_1 = find_index(current_z, &z_vals_1);
    index_2 = find_index(current_z, &z_vals_2);
    rebinned_data[1][z] = kernel_1[index_1];
    rebinned_data[2][z] = kernel_2[index_2];
    
    if(current_z >= z_vals_1[Nz_1-1] || current_z < z_vals_1[0]) rebinned_data[1][z] = 0.0;
    if(current_z >= z_vals_2[Nz_2-1] || current_z < z_vals_2[0]) rebinned_data[2][z] = 0.0;
  }
  
  return rebinned_data;
  
}
