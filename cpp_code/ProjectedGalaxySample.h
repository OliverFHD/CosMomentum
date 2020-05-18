#include <gsl/gsl_sf_gamma.h>


class ProjectedGalaxySample : public GalaxySample {

 public:

   ProjectedGalaxySample(Matter* matter, double density_in_arcmin_squared, double b1, double b2, double a0, double a1, string input_file);
   ~ProjectedGalaxySample();
   
   void set_n_of_z_data(string input_file);
   void set_parameters_projected(double density_in_arcmin_squared, double b1, double b2, double a0, double a1, string input_file);
   double set_b2_to_minimise_negative_densities_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale);
   double compute_variance_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale);
   
   int return_N_max_in_angular_tophat(double theta_in_arcmin, double var_NL_rescale);
   vector<double> return_CiC_PDF_in_angular_tophat(double theta_in_arcmin, double f_NL, double var_NL_rescale);
   
 private:
   
   vector<double> z_values;
   vector<double> n_of_z_values;
   
   double density; // in radians
   
};
