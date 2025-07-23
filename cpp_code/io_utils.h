
using namespace std;

vector<vector<double> > read_table_double(string file_name){
  
  string line;
  double dummy;
  
  fstream in;
  in.open(file_name);
  
  vector<vector<double> > matrix(0,vector<double>(0,0.0));
  
  while(in.good()){
    getline(in, line);
    if(line[0]!='#' && in.good()){
      vector<double> row(0, 0.0);
      stringstream sstr;
      sstr << line;
      while(sstr.good()){
        sstr >> dummy;
        row.push_back(dummy);
      }
      matrix.push_back(row);
    }
  }
  
  return matrix;
  
}



vector<vector<double> > read_table_double_transposed(string file_name){
  
  string line;
  double dummy;
  
  fstream in;
  in.open(file_name);
  
  vector<vector<double> > matrix(0,vector<double>(0,0.0));
  
  while(in.good()){
    getline(in, line);
    if(line[0]!='#' && in.good()){
      vector<double> row(0, 0.0);
      stringstream sstr;
      sstr << line;
      while(sstr.good()){
        sstr >> dummy;
        row.push_back(dummy);
      }
      matrix.push_back(row);
    }
  }
  
  int N_row = matrix.size();
  int N_col = matrix[0].size();
  
  vector<vector<double> > matrix_transposed(N_col,vector<double>(N_row,0.0));
  
  for(int i = 0; i < N_col; i++){
    for(int j = 0; j < N_row; j++){
      matrix_transposed[i][j] = matrix[j][i];
    }
  }
  
  return matrix_transposed;
  
}
