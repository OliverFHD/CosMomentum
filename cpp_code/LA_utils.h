
using namespace std;

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
