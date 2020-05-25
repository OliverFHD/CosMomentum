#ifndef _FriedmannUniverse
#define _FriedmannUniverse

using namespace constants;



class FriedmannUniverse {

 public:

  FriedmannUniverse();
  ~FriedmannUniverse();
  
  int return_number_of_time_steps();
  double return_a_initial();
  double return_a_final();
  double return_eta_initial();
  double return_eta_final();
  
  void print_background_cosmology(string filename);
  
  double f_k(double w);
  double t_at_eta(double e);
  double a_at_eta(double e);
  double H_at_eta(double e);
  double H_prime_at_eta(double e);  
  double eta_at_a(double a);  
  double eta_i(int i);
  
  void set_spatial_curvature(double curv){this->spatial_curvature_today = curv;};
  
  vector<vector<double> > return_background_expansion();
  vector<vector<double> > return_background_expansion(int conformal_time_steps);

 private:
   
   double spatial_curvature_today; // in units with c/H_0 = 1

 protected:
   
   vector<double> t;
   vector<double> eta;
   vector<double> a;
   vector<double> H;
   vector<double> H_prime;
   
};

#include "FriedmannUniverse.cpp"
#endif
