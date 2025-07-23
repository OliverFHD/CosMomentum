#include <vector>
#include <gsl/gsl_sf_bessel.h>


using namespace std;
using namespace constants;


struct cosmological_model{

  int collapse;
  
  double Omega_m;
  double Omega_r;
  double Omega_L;
  double Omega_b;
  double Omega_k;
  
  double theta_27;
  double w0;
  double w1;
  double h_100;
  
  double n_s;
  double sigma_8;

};



void give_area_from_cosines(double cos_th_min, double cos_th_max, double* A){
  (*A) = constants::pi2*(cos_th_min-cos_th_max);
}



double w_1(double k){
  
  return 3.0*(sin(k)-k*cos(k))/pow(k,3);
  
}

double deriv_of_w_1(double k){
  
  return 9.0*(k*cos(k)+(k*k/3.0-1.0)*sin(k))/pow(k,4);
  
}

double w_R(double k, double R){
  
  return w_1(k*R);
  
}

double deriv_of_w_R(double k, double R){
  
  return k*deriv_of_w_1(k*R);
  
}

double w_R_2D(double k, double R){
  
  double x = k*R;
	if(x > eps_from_0)
		return 2.0*gsl_sf_bessel_J1(x)/x;

  return 1.0;
}

double w_L_1D(double k, double L){
  
  double x = k*L;
  if(x > eps_from_0)
    return 2.0*sin(x/2.0)/x;

  return 1.0;
}

double w_L_1D_Gauss(double k, double L){
  
  double x = k*L;
  if(x > eps_from_0)
    return exp(-0.5*pow(x/sqrtPi_times_2,2));

  return 1.0;
}

double deriv_of_w_R_2D(double k, double R){
  
  double x = k*R;
  return -2.0*gsl_sf_bessel_Jn(2, x)/x*k;
  
}

double second_deriv_of_w_R_2D(double k, double R){
  
  double x = k*R;
  double one_over_x = 1.0/x;
  return k*k*(2.0*gsl_sf_bessel_Jn(2, x)*pow(one_over_x,2) + one_over_x*gsl_sf_bessel_Jn(3, x) - one_over_x*gsl_sf_bessel_J1(x));
  
}
















/*******************************************************************************************************************************************************
 * These functions are for the Smith_et_al fitting formula
*******************************************************************************************************************************************************/




double a_n(double n_spec, double C, double Om_w, double w){
 
  //double x = 1.4861 + (1.8369 + (1.6762 + (0.7940 + 0.1670*n_spec)*n_spec)*n_spec)*n_spec - 0.6206*C;
  double x = 1.5222 + (2.8553 + (2.3706 + (0.9903 + 0.2250*n_spec)*n_spec)*n_spec)*n_spec - 0.6038*C + 0.1749*Om_w*(1.0 + w);
  
  return pow(10.0, x);
  
}

double b_n(double n_spec, double C, double Om_w, double w){
 
  //double x = 0.9463 + (0.9466 + 0.3084*n_spec)*n_spec  -0.94*C;
  double x = -0.5642 + (0.5864 + 0.5716*n_spec)*n_spec  - 1.5474*C + 0.2279*Om_w*(1.0 + w);
  
  return pow(10.0, x);
  
}

double c_n(double n_spec, double C){
 
  //double x = -0.2807 + (0.6669 + 0.3214*n_spec)*n_spec -0.0793*C;
  double x = 0.3698 + (2.0404 + 0.8161*n_spec)*n_spec + 0.5869*C;
  
  return pow(10.0, x);
  
}

double mu_n(double n_spec){
 
  //double x = -3.5442 + 0.1908*n_spec;  
  //return pow(10.0, x);
  return 0.0;
  
}

double nu_n(double n_spec){

  //double x = 0.9589 + 1.2857*n_spec;
  double x = 5.2105 + 3.6902*n_spec;
  
  return pow(10.0, x);
  
}

double gamma_n(double n_spec, double C){
 
  //return 0.8649 + 0.2989*n_spec + 0.1631*C;
  return 0.1971 - 0.0843*n_spec + 0.8460*C;
  
}

//double alpha_n(double n_spec){
double alpha_n(double n_spec, double C){
 
  //return 1.3884 + (0.3700 - 0.1452*n_spec)*n_spec;
  return abs(6.0835 + (1.3373 - 0.1959*n_spec)*n_spec - 5.5274*C);
  
}

//double beta_n(double n_spec){
double beta_n(double n_spec, double C){
 
  //return 0.8291 + (0.9854 + 0.3401*n_spec)*n_spec;
  return 2.0379 + ( -0.7354 + (0.3157 + (1.2490 + 0.3980*n_spec)*n_spec)*n_spec)*n_spec - 0.1682*C;
  
}

double f_1a(double om_m){
  return pow(om_m, -0.0732);
}

double f_2a(double om_m){
  return pow(om_m, -0.1423);
}

double f_3a(double om_m){
  return pow(om_m, 0.0725);
}


double f_1b(double om_m){
  return pow(om_m, -0.0307);
}

double f_2b(double om_m){
  return pow(om_m, -0.0585);
}

double f_3b(double om_m){
  return pow(om_m, 0.0743);
}

double f_no_index(double x){
  return x/4.0*(1+x/2.0);
}





vector<double> factoria(int n){
 
  vector<double> fact(1, 1.0);
  
  for(int i = 1; i <= n; i++){
    fact.push_back(fact[i-1]*double(i));
  }
  
  return fact;
  
}


extern "C" double Bell_polynomial(int n, int k, vector<double> *x, vector<double> *fact){
   
  if(n == 0 && k == 0) return 1.0;
  if(n == 0 && k != 0) return 0.0;
  if(n != 0 && k == 0) return 0.0;
  
  
  double sum = 0.0;
  double binomial;
  
  
  for(int i = 1; i <= n-k+1; i++){
    binomial = (*fact)[n-1]/(*fact)[n-i]/(*fact)[i-1];
    sum += (*x)[i]*binomial*Bell_polynomial(n-i, k-1, x, fact);
  }
  
  
  return sum;
  
}


void Faa_di_Bruno(vector<double> *Gy_coeffs, vector<double> *Gtau_coeffs, vector<double> *y_coeffs){
  
  int n = Gtau_coeffs->size();
  
  vector<double> fact = factoria(n);
  vector<double> x(n, 0.0);
  
  vector<vector<double> > Bell_matrix(n, vector<double>(n, 0.0));
  vector<vector<double> > inverse_Bell_matrix(n, vector<double>(n, 0.0));
  
  for(int i = 0; i < n; i++){
    x[i] = (*y_coeffs)[i]*fact[i];
  }
  
  
  
  /*
   * Gtau_coeffs[i] = sum_j Bell_matrix[i][j] * Gy_coeffs[i];
   * 
   * => need inverse matrix later
   */
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j <= i; j++){
      Bell_matrix[i][j] = Bell_polynomial(i, j, &x, &fact);
    }
  }
  
  invert_matrix(&Bell_matrix, &inverse_Bell_matrix);
  
  (*Gy_coeffs) = vector<double>(n, 0.0);
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j <= i; j++){
      (*Gy_coeffs)[i] += inverse_Bell_matrix[i][j]*(*Gtau_coeffs)[j]*fact[j];
    }
    (*Gy_coeffs)[i] /= fact[i];
  }
  
}


void Faa_di_Bruno_inverse(vector<double> *Gtau_coeffs, vector<double> *Gy_coeffs, vector<double> *y_coeffs){
  
  int n = Gy_coeffs->size();
  
  vector<double> fact = factoria(n);
  vector<double> x(n, 0.0);
  
  vector<vector<double> > Bell_matrix(n, vector<double>(n, 0.0));
  
  for(int i = 0; i < n; i++){
    x[i] = (*y_coeffs)[i]*fact[i];
  }

  
  for(int i = 0; i < n; i++){
    for(int j = 0; j <= i; j++){
      Bell_matrix[i][j] = Bell_polynomial(i, j, &x, &fact);
    }
  }
  
  
  (*Gtau_coeffs) = vector<double>(n, 0.0);
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j <= i; j++){
      (*Gtau_coeffs)[i] += Bell_matrix[i][j]*(*Gy_coeffs)[j]*fact[j];
    }
    (*Gtau_coeffs)[i] /= fact[i];
  }
  
}


void return_Bell_matrix(vector<vector<double> > *Bell_matrix, vector<vector<double> > *inverse_Bell_matrix, vector<double> *y_coeffs){
  
  int n = y_coeffs->size();
  
  vector<double> fact = factoria(n);
  vector<double> x(n, 0.0);
  
  Bell_matrix->resize(n, vector<double>(n, 0.0));
  inverse_Bell_matrix->resize(n, vector<double>(n, 0.0));
  
  for(int i = 0; i < n; i++){
    x[i] = (*y_coeffs)[i]*fact[i];
  }
  
  
  
  /*
   * Gtau_coeffs[i] = sum_j Bell_matrix[i][j] * Gy_coeffs[i];
   * 
   * => need inverse matrix later
   */
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j <= i; j++){
      (*Bell_matrix)[i][j] = Bell_polynomial(i, j, &x, &fact);
    }
  }
  
  invert_matrix(Bell_matrix, inverse_Bell_matrix);
  
}









