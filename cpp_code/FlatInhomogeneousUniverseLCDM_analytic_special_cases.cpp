
using namespace std;

/*
 * The Friedmann equations in terms of conformal quantities are:
 * 
 *      H′ = -4pi/3*(rho + 3p)*a^2
 * H^2 + k = 8pi/3*rho*a^2
 * 
 * Or, measuring H in units of H_0 and eta in units of 1/H_0 and defining Omega_k = -k:
 * 
 * 
 *      H′ = -0.5*(Omega_m/a + 2*Omega_r/a^2 - 2*Omega_L*a^2)
 *     H^2 = (Omega_m/a + Omega_r/a^2 + Omega_k + Omega_L*a^2)
 * 
 * NOTE: This module is outputting the growing mode of linear structure growth only.
 * 
 * NOTE: Technically, there must be a non-zero fraction of matter density in order to have density fluctuations. This means that there will always be redshifts, at which the following calculations are incorrect.
 * 
 */


/*
 * In an EdS universe in units of H_0 we have:
 * 
 *    H_conformal   = 2.0/eta
 * 
 *    H_conformal^2 = 1.0/a
 * 
 * =>             a = 1.0/H_conformal^2
 *                  = eta^2/4.0
 * 
 *     H_conformal' = -0.5/a
 *                  = -0.5*H_conformal^2
 *                  = -2.0/eta^2
 * 
 *          t_phys' = a
 *                  = eta^2/4.0
 * =>        t_phys = eta^3/12.0
 * 
 */
void FlatInhomogeneousUniverseLCDM::growth_of_DM_fluctuations_in_flat_matter_dominated_universe(double a, double *eta, double *D, double *D_prime){
  {
    using namespace error_handling;
    string error_message, warning_message;
    
    error_message = "Expansion rate ill defined at a <= 0.0 .";
    test_value_of_double(a, 0.0, LARGER, error_message);
  }
  
  double H_conformal = 1.0/sqrt(a);
  
  
  (*eta) = 2.0*sqrt(a);
  (*D) = a;
  (*D_prime) = H_conformal*a;
  
}


/*
 * In a flat, radiation dominated universe in units of H_0 we have:
 * 
 *    H_conformal   = 1.0/eta
 * 
 *    H_conformal^2 = 1.0/a^2
 * 
 * =>             a = 1.0/H_conformal
 *                  = eta
 * 
 *     H_conformal' = -1/a^2
 *                  = -1/eta^2
 * 
 *          t_phys' = a
 *                  = eta
 * =>        t_phys = eta^2/2.0
 * 
 */

void FlatInhomogeneousUniverseLCDM::growth_of_DM_fluctuations_in_flat_radiation_dominated_universe(double a, double *eta, double *D, double *D_prime){
  {
    using namespace error_handling;
    string error_message, warning_message;
    
    error_message = "Expansion rate ill defined at a <= 0.0 .";
    test_value_of_double(a, 0.0, LARGER, error_message);
  }
  
  (*eta) = a;
  (*D) = 1.0;
  (*D_prime) = 0.0;
  
}


/*
 * In a flat de-Sitter universe in units of H_0 we have:
 * 
 *    H_conformal^2 = a^2
 * 
 * =>             a = -1.0/eta
 *      H_conformal = -1.0/eta
 *     H_conformal' = 1/eta^2
 * 
 *          t_phys' = a
 *                  = -1.0/eta
 * =>        t_phys = -ln(-eta) , since eta < 0
 * 
 */
void FlatInhomogeneousUniverseLCDM::growth_of_DM_fluctuations_in_flat_Lambda_dominated_universe(double a, double *eta, double *D, double *D_prime){
  {
    using namespace error_handling;
    string error_message, warning_message;
    
    error_message = "Expansion rate ill defined at a <= 0.0 .";
    test_value_of_double(a, 0.0, LARGER, error_message);
  }
  
  (*eta) = -1.0/a;
  (*D) = pow((*eta), 2.0)/2.0;
  (*D_prime) = (*eta);
}
