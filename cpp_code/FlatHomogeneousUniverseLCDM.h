#ifndef _FlatHomogeneousUniverseLCDM
#define _FlatHomogeneousUniverseLCDM
#include "FriedmannUniverse.h"


using namespace constants;


class FlatHomogeneousUniverseLCDM : public FriedmannUniverse {

 public:

  FlatHomogeneousUniverseLCDM(cosmological_model cosmo, double a_min, double a_max);
  ~FlatHomogeneousUniverseLCDM();
   
  double rho_m_of_a(double scale); // All in units of TODAYS critical density
  double rho_r_of_a(double scale);
  double rho_L_of_a(double scale);
  double w_L_of_a(double scale);
  
  double return_Omega_m(){return this->cosmology.Omega_m;};
  double return_Omega_b(){return this->cosmology.Omega_b;};
  double return_Omega_r(){return this->cosmology.Omega_r;};
  double return_Omega_L(){return this->cosmology.Omega_L;};
  double return_h_100(){return this->cosmology.h_100;};
  double return_theta_27(){return this->cosmology.theta_27;};
  cosmological_model return_cosmology();
  
  double Mass_within_R_in_Mpc_over_h(double R_in_Mpc_over_h);
  double R_in_Mpc_over_h_enclosing_M(double M_in_Msol_over_h);
  
  static void expansion_in_flat_matter_dominated_universe(double a, double Omega_m, double *t_phys, double *eta, double *H_conformal, double *H_conformal_prime);
  static void expansion_in_flat_radiation_dominated_universe(double a, double *t_phys, double *eta, double *H_conformal, double *H_conformal_prime);
  static void expansion_in_flat_Lambda_dominated_universe(double a, double *t_phys, double *eta, double *H_conformal, double *H_conformal_prime);
  
 protected:
   
   double t_initial;
   double eta_initial;
   double H_initial;
   double H_prime_initial;
   double a_initial;
   double a_final;
   
   cosmological_model cosmology;
   
 private:
   
   void set_initial_conditions();
   void set_expansion_history();
   
   double hubble_from_Friedmann(double a_start);
   double hubble_prime_from_Friedmann(double a_start);
   
};

#include "FlatHomogeneousUniverseLCDM.cpp"
#endif
