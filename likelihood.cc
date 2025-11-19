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
  double resolution = 0.1;
  double starting_value = 0;
  double stoping_value = 6;
  double unc_method1;
  
  
  
  raw_data.open("datensumme.txt");
  output.open("likelihood.txt");
  output2.open("nll.txt");
  output3.open("deltanll.txt");
  
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
  //cout << "statistic uncertainty:"<<unc_method1<<endl;
  //cout<< " calculated uncertainty:"<<new_count<<endl;
  cout << "likelihood: "<<likelihood<<endl;

  
  raw_data.close();
  output.close();
  output2.close();
  output3.close();


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
  pois;
  return pois;
}
 
