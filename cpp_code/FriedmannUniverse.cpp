#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include "cosmology_utils.cpp"
 
using namespace std;
using namespace constants;


/*******************************************************************************************************************************************************
 * FriedmannUniverse::FriedmannUniverse
 * 
 * Constructor of class FriedmannUniverse
 * 
*******************************************************************************************************************************************************/


FriedmannUniverse::FriedmannUniverse(){
  
}

/*******************************************************************************************************************************************************
 * FriedmannUniverse::~FriedmannUniverse
 * 
 * Destructor of class FriedmannUniverse
 * 
 * 
*******************************************************************************************************************************************************/

FriedmannUniverse::~FriedmannUniverse(){
  
}


int FriedmannUniverse::return_number_of_time_steps(){
  return this->eta.size();
}

/*******************************************************************************************************************************************************
 * FriedmannUniverse::return_eta_initial
 * 
 * Return lowest value of eta for which expansion history is available.
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::return_eta_initial(){
    
  return this->eta[0];
  
}

/*******************************************************************************************************************************************************
 * FriedmannUniverse::eta_i
 * 
 * Return ith value of eta for which expansion history is available.
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::eta_i(int i){
  if(i > this->return_number_of_time_steps()-1){
    using namespace error_handling;
    general_error_message("i = " + to_string(i) + " should be < " + to_string(this->return_number_of_time_steps()) + " (number of time steps stored in this->eta) in method FriedmannUniverse::eta_i(int i).");
  }
  return this->eta[i];
}



/*******************************************************************************************************************************************************
 * FriedmannUniverse::return_a_initial
 * 
 * Return lowest value of scale factor for which expansion history is available.
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::return_a_initial(){
    
  return this->a[0];
  
}


/*******************************************************************************************************************************************************
 * FriedmannUniverse::return_eta_final
 * 
 * Return largest value of eta for which expansion history is available.
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::return_eta_final(){
    
  return this->eta[this->eta.size()-1];
  
}

/*******************************************************************************************************************************************************
 * FriedmannUniverse::return_a_final
 * 
 * Return largest value of scale factor for which expansion history is available.
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::return_a_final(){
    
  return this->a[this->a.size()-1];
  
}

/*******************************************************************************************************************************************************
 * f_k
 * Description:
 *  - Co-moving angular diameter distance.
 * Arguments:
 *  - w: co-moving distance at which f_k shall be evaluated
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::f_k(double w){
  
  
  if(abs(this->spatial_curvature_today) < eps_from_0)
    return w;
  
  double curvature_radius_today = 1.0/sqrt(abs(this->spatial_curvature_today));

  if(this->spatial_curvature_today > 0.0)
    return curvature_radius_today*sin(w/curvature_radius_today);
  
  return curvature_radius_today*sinh(w/curvature_radius_today);

  
}

/*******************************************************************************************************************************************************
 * 3.9 a_at_eta
 * Description:
 *  - Gives value of scale factor at given conformal time.
 * Arguments:
 *  - e: conformal time
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::a_at_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->a, constants::order_of_interpolation);
  
}



/*******************************************************************************************************************************************************
 * 3.9 a_at_eta
 * Description:
 *  - Gives value of scale factor at given conformal time.
 * Arguments:
 *  - e: conformal time
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::t_at_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->t, constants::order_of_interpolation);
  
}

/*******************************************************************************************************************************************************
 * 3.10 H_at_eta
 * Description:
 *  - Gives value of conformal expansion rate at given conformal time.
 * Arguments:
 *  - e: conformal time
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::H_at_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->H, constants::order_of_interpolation);
  
}

/*******************************************************************************************************************************************************
 * 3.11 H_prime_at_eta
 * Description:
 *  - Gives derivative of conformal expansion rate at given conformal time.
 * Arguments:
 *  - e: conformal time
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::H_prime_at_eta(double e){
  
  return interpolate_neville_aitken(e, &this->eta, &this->H_prime, constants::order_of_interpolation);
  
}

/*******************************************************************************************************************************************************
 * 3.12 eta_at_a
 * Description:
 *  - Gives value of conformal time at given scale factor.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/

double FriedmannUniverse::eta_at_a(double a){
  
  return interpolate_neville_aitken(a, &this->a, &this->eta, constants::order_of_interpolation);
  
}


/*******************************************************************************************************************************************************
 * 3.13 eta_at_a
 * Description:
 *  - Gives value of conformal time at given scale factor.
 * Arguments:
 *  - a: scale factor
 * 
*******************************************************************************************************************************************************/


vector<vector<double> > FriedmannUniverse::return_background_expansion(){
  
  int n_column = 5;
  int time_steps = this->a.size();
  vector<vector<double> > background_expansion(n_column, vector<double>(time_steps, 0.0));
  
  background_expansion[0] = this->t;
  background_expansion[1] = this->eta;
  background_expansion[2] = this->a;
  background_expansion[3] = this->H;
  background_expansion[4] = this->H_prime;
  
  return background_expansion;
  
}

vector<vector<double> > FriedmannUniverse::return_background_expansion(int conformal_time_steps){
  
  int n_column = 5;
  vector<vector<double> > background_expansion(n_column, vector<double>(conformal_time_steps, 0.0));
  
  double a_max = this->a[this->a.size()-1];
  double a_min = this->a[0];
  double da = (a_max - a_min)/double(conformal_time_steps-1);
  double a;
  double e;
  
  for(int i = 0; i < conformal_time_steps; i++){
    a = a_min + double(i)*da;
    e = eta_at_a(a);
    background_expansion[0][i] = this->t_at_eta(e);
    background_expansion[1][i] = e;
    background_expansion[2][i] = a;
    background_expansion[3][i] = this->H_at_eta(e);
    background_expansion[4][i] = this->H_prime_at_eta(e);
  }
  
  return background_expansion;
  
}



/*******************************************************************************************************************************************************
 * 4.1 print_background_cosmology
 * Description:
 *  - Prints table with expansion history of the universe.
 * Arguments:
 *  - file_name: file name to which the table shall be output
 * 
*******************************************************************************************************************************************************/

void FriedmannUniverse::print_background_cosmology(string file_name){
  
  fstream cosmo_stream;
  
  remove(file_name.c_str());
  FILE * F = fopen(file_name.c_str(), "w");
  fclose(F);
  cosmo_stream.open(file_name.c_str());
  
  cosmo_stream << setprecision(10) << scientific;
  cosmo_stream << "#a_max = " << this->return_a_final() << '\n';
  cosmo_stream << "#t" << setw(20) << "eta(t)" << setw(20) << "w(eta)" << setw(20) << "z(eta)" << setw(20) << "a(eta)" << setw(20) << "H(eta)" << setw(20) << "H_prime(eta)\n";
  for(int i = 0; i < this->a.size(); i++){
    cosmo_stream << this->t[i];
    cosmo_stream << setw(20) << this->eta[i];
    cosmo_stream << setw(20) << this->eta_at_a(1.0)-this->eta[i];
    cosmo_stream << setw(20) << 1.0/this->a[i] - 1.0;
    cosmo_stream << setw(20) << this->a[i];
    cosmo_stream << setw(20) << this->H[i];
    cosmo_stream << setw(20) << this->H_prime[i] << '\n';
  }
  
  cosmo_stream.close();
  
}




  




