//moin!

#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<numeric>


using namespace std;


double poisson(double, int);
double mean_calculator(vector<int>);
double prob(vector<int>, double);
double uncertainty_calculation(vector<int>);
long double lambda_calculation(vector<int>, double);


int main() {
  ifstream raw_data;
  ofstream output;
  ofstream output2, output3;
  raw_data.open("datensumme.txt");
  output.open("likelihood.txt");
  output2.open("nll.txt");
  output3.open("deltanll.txt");
  vector<int> total_data;
  int col1;
  double starting_value = 0;
  double stoping_value = 6;
  //4.1a)
  while(raw_data>>col1)
    total_data.push_back(col1);
  double mean;
  mean = mean_calculator(total_data);
  double likelihood;
  likelihood = prob(total_data, mean);
  double best_log_like;
  best_log_like = -2 * log(likelihood);
  
  //4.1b,c,d)
  double est_mean;
  est_mean = starting_value;
  double result;
  double ln_results;
  double delta;
  int count = 0;
  double resolution = 0.1;
  while(est_mean < stoping_value){
    result = prob(total_data, est_mean);
    ln_results = -2 * log(result);
    delta = ln_results - best_log_like;
    if (delta < 1.0)
      count++;
    vector<double> likelihoods;
    likelihoods.push_back(result);
    output<<est_mean<<'\t'<<result<<endl;
    output2<<est_mean<<'\t'<<ln_results<<endl;
    output3<<est_mean<< '\t'<<delta<<endl;
    est_mean +=  resolution;
  }
  double new_count;
  new_count = double(count) * (resolution);
  double unc_method1;
  unc_method1  = uncertainty_calculation(total_data);
  //cout << "statistic uncertainty:"<<unc_method1<<endl;
  //cout<< " calculated uncertainty:"<<new_count<<endl;
  cout << "likelihood: "<<likelihood<<endl;

  //4.1e
  long double lambda;
  lambda = lambda_calculation(total_data, mean);
  cout <<"calculated lambda: "<< lambda<<endl;
  cout <<"calculated -2 ln lambda: "<< -2*log(lambda)<<endl;
  double ndof = 233;
  double z;
  z = (-2*log(lambda)-ndof)/sqrt(2*ndof);
  cout<<"calculated z: "<< z<<endl;

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

  
double prob(vector<int> zahlen, double mean){
  double likelihood;
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

long double lambda_calculation(vector<int> zahlen, double mean){
    long double lambda;
    lambda = 1;
    for (int k : zahlen){
        lambda *= poisson(mean, k)/poisson(k, k);
  }
  return lambda;
}
 
