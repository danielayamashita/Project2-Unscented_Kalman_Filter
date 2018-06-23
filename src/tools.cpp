#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
   VectorXd res;
  VectorXd rmse(4);
  VectorXd sum(4);
  rmse << 0,0,0,0;
  sum << 0,0,0,0;

   // check the validity of the following inputs:
  if ((estimations.size()==0) || (estimations.size() != ground_truth.size()))
  {
	    cout<< "Error: estimation or ground_truth does not have a good size."<< endl;
	    return rmse;
	}
  //accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code 
		res = (estimations[i] - ground_truth[i]);
		res = res.array()*res.array();
		rmse = rmse+res; 
	}
  rmse/=estimations.size();
	rmse = rmse.array().sqrt();

  //std::cout << "RMSE:" <<rmse<< std::endl;
  return rmse;
}