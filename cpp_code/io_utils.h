
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
