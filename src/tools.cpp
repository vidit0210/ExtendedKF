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

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;


  //The estimation vector size should not be zero
  if (estimations.size() == 0) {
    cout << "Estimation size is zero" << endl;
    return rmse;
  }
  else
  {
    //The estimation vector size should be equal to ground truth vector size
    if (estimations.size() != ground_truth.size())
    {
      cout << "Invalid ground_truth data" << endl;
      return rmse;
    }
  }

  for (unsigned int i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // return squared root of mean RMSE
  return (rmse / estimations.size()).array().sqrt();;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);

  //Recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);


  //Pre-compute a set of terms to avoid repeated calculation
  float c1 = px * px + py * py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //Check division by zero
  if (fabs(c1) < 0.0001) {
    std::cout << "ERROR:: CalculateJacobian () -- Division by Zero" << std::endl;
    Hj.fill(0.0);
    return Hj;
  }

  //Compute the Jacobian matrix
  Hj << (px / c2), (py / c2), 0, 0,
    -(py / c1), (px / c1), 0, 0,
    py*(vx*py - vy * px) / c3, px*(px*vy - py * vx) / c3, px / c2, py / c2;

  return Hj;
}
