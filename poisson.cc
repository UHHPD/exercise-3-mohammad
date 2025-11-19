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
  vector<int> total_data;
  vector<long double> likelihoods;
  ifstream raw_data;
  ofstream output;
  ofstream output2, output3;
  int col1;
  double mean;
  long double likelihood;
  double best_log_like;
  double est_mean;
  long double result;
  long double ln_results;
  long double delta;
  int count = 0;
  double new_count;
  double resolution = 0.001;
  double starting_value = 1;
  double stoping_value = 6;
  double unc_method1;
  vector<int> zaehler(11);
  vector<double> poisson_estimations(11);
  ifstream summs_file;
  ofstream previous_output_file;
  ofstream previous_output_file2;

  
  
  raw_data.open("datensumme.txt");
  output.open("likelihood.txt");
  output2.open("nll.txt");
  output3.open("deltanll.txt");
  previous_output_file.open("hist.txt");
  previous_output_file2.open("histpoi.txt");


  zaehler = abundance_calculator(raw_data);
  poisson_estimations = poisson_collector(zaehler);

  raw_data.clear();
  raw_data.seekg(0, ios::beg);

  
  for(int i=0; i< zaehler.size(); i++){
    //cout<< i << " : " << zaehler[i]<<endl;
    previous_output_file << i << '\t' << zaehler[i] << endl;
    previous_output_file2 << i << '\t' << zaehler[i] << '\t' << poisson_estimations[i] << endl;
    cout << i << '\t' << zaehler[i] <<endl; //'\t' << poisson_estimations[i] << endl;
  }


  
  while(raw_data>>col1)
    total_data.push_back(col1);
  
  mean = mean_calculator(total_data);
  likelihood = prob(total_data, mean);
  best_log_like = -2 * log(likelihood);
  

  est_mean = starting_value;
  while(est_mean < stoping_value){
    result = prob(total_data, est_mean);
    ln_results = -2 * log(result);
    delta = ln_results - best_log_like;
    if (delta < 1.0)
      count++;
    likelihoods.push_back(result);
    output<<est_mean<<'\t'<<result<<endl;
    output2<<est_mean<<'\t'<<ln_results<<endl;
    output3<<est_mean<< '\t'<<delta<<endl;
    est_mean +=  resolution;
  }
  new_count = double(count) * (resolution);
  unc_method1  = uncertainty_calculation(total_data);
  cout << "statistic uncertainty:"<<unc_method1<<endl;
  cout<< " calculated uncertainty:"<<new_count<<endl;


  
  raw_data.close();
  output.close();
  output2.close();
  output3.close();
  previous_output_file.close();
  previous_output_file2.close();


  return 0;
}


double uncertainty_calculation(vector<int> zahlen){
  double mean;
  double unc;
  mean = mean_calculator(zahlen);
  unc = mean / sqrt(zahlen.size());
  return unc;
}

double mean_calculator(vector<int> zahlen){
      double sum = accumulate(zahlen.begin(), zahlen.end(), 0.0);
      double mean = sum / zahlen.size();
      return mean;
}

  
long double prob(vector<int> zahlen, double mean){
  long double likelihood;
  double poisson_value;
  likelihood = 1;
  for (int k : zahlen){
    poisson_value = poisson(mean, k);
    likelihood *= poisson_value;
  }
  return likelihood;
}
    
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
