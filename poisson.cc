//moin!

#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<numeric>


using namespace std;


double poisson(double, int);
double mean_calculator(vector<int>);
long double prob(vector<int>, double);
double uncertainty_calculation(vector<int>);
vector<int> abundance_calculator(ifstream&);
vector<double> poisson_collector(vector<int>);

int main() {
  ifstream raw_data;
  raw_data.open("datensumme.txt");
  ofstream output,output2, output3;
  output.open("likelihood.txt");
  output2.open("nll.txt");
  output3.open("deltanll.txt");
  ofstream previous_output_file;
  previous_output_file.open("hist.txt");
  ofstream previous_output_file2;
  previous_output_file2.open("histpoi.txt");

  vector<int> zaehler(11);
  zaehler = abundance_calculator(raw_data);
  vector<double> poisson_estimations(11);
  poisson_estimations = poisson_collector(zaehler);

  raw_data.clear();
  raw_data.seekg(0, ios::beg);

  // 3.1abc)
  for(int i=0; i< zaehler.size(); i++){
    //cout<< i << " : " << zaehler[i]<<endl;
    previous_output_file << i << '\t' << zaehler[i] << endl;
    previous_output_file2 << i << '\t' << zaehler[i] << '\t' << poisson_estimations[i] << endl;
    cout << i << '\t' << zaehler[i] <<endl; //'\t' << poisson_estimations[i] << endl;
  }

  return 0;
}

    
// 3.1c)
double poisson(double mean, int observation){
  int gamma;
  double mu_ka, exp_mu, pois;
  gamma = tgamma(double(observation + 1));
  mu_ka = pow(mean, observation);
  exp_mu = exp(-1 * mean);
  pois = (mu_ka * exp_mu) / (gamma);
  return pois;
}

vector<int> abundance_calculator(ifstream& file){
  vector<int> abund_vect(11);
  int col1;

  while(file >> col1)
    abund_vect[col1] += 1;
  return abund_vect;
}


vector<double> poisson_collector(vector<int> hists){
  vector<double> poisson_values(hists.size());
  double mean, sum_games, sum_goals;

  sum_goals = 0.;
  sum_games = accumulate(hists.begin(), hists.end(), 0.0);
  for(int i=0; i < hists.size(); i++)
    sum_goals += (hists[i] * i);
  
  mean = sum_goals / sum_games;  
  
  for(int i=0; i < hists.size(); i++)
    poisson_values[i] = poisson(mean, i) * sum_games;

  return poisson_values;
}
