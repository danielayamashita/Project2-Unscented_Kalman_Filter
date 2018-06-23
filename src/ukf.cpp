#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //set state dimension
  n_x = 5;

  //set augmented dimension
  n_aug = 7;

  Xsig_pred_ = MatrixXd(n_x, 2 * n_aug + 1);

  time_us_ = 0;

  is_initialized_ = false;

  // Vector for weights
  weights_ = VectorXd(15);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) 
  {
    #ifdef PRINT_HEADERS_H_
      std::cout << "--------------------------------------" << std::endl;
      std::cout << "Process Measurement" << std::endl;

    #endif
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      #ifdef PRINT_HEADERS_H_
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "INIT RADAR Measurement:" << std::endl;
        std::cout << "[*]:" <<meas_package.raw_measurements_<< std::endl;
        std::cout << "[0]: "<<meas_package.raw_measurements_[0] << std::endl;
        std::cout << "[0]: "<<meas_package.raw_measurements_[1] << std::endl;
      #endif
      x_(0) = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
      x_(1) = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);

      x_(2) =  meas_package.raw_measurements_[2]*cos(meas_package.raw_measurements_[1]);
      x_(3) =  meas_package.raw_measurements_[2]*sin(meas_package.raw_measurements_[1]);
      x_(5) = 0;
      if(x_(0) <0.00001)
      {
        x_(0) = 0.00001;
      }
      if(x_(1) <0.00001)
      {
        x_(1) = 0.00001;
      }
      // Inicialize the covariance matrix P
      //state covariance matrix
      P_ << std_radr_*std_radr_, 0, 0, 0, 0,
            0, std_radr_*std_radr_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, std_radphi_, 0,
            0, 0, 0, 0, std_radphi_;
      
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) 
    {
      #ifdef PRINT_HEADERS_H_
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "INIT LASER Measurement:" << std::endl;
        std::cout << "[*]:" <<meas_package.raw_measurements_<< std::endl;
        std::cout << "[0]:" <<meas_package.raw_measurements_[0]<< std::endl;
        std::cout << "[1]:" <<meas_package.raw_measurements_[1]<< std::endl;
      #endif
      x_<< meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0,0;
     #ifdef PRINT_HEADERS_H_
        std::cout << "x" <<x_<< std::endl;
      #endif
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    }




    #ifdef PRINT_HEADERS_H_
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "End ProcessMeasurement" << std::endl;
    std::cout << "timestamp:" <<time_us_<< std::endl<< std::endl;
    std::cout << "P:" <<P_<< std::endl;
    #endif 
    time_us_ = meas_package.timestamp_ ;
    
    is_initialized_ = true;
    return;
 }

 
 double dt = (meas_package.timestamp_ - time_us_)/1000000.0;	//dt - expressed in seconds
 time_us_ = meas_package.timestamp_ ;


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) 
  {
    // Radar update
    UpdateRadar(meas_package);
  } 
  else 
  {
    // Laser updates
    UpdateLidar(meas_package);
  }


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //---------------------------------------------------------------------------------
  //Generate Sigma Points
  //---------------------------------------------------------------------------------
  MatrixXd Xsig = MatrixXd(5, 11);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  lambda_ = 3 - n_x_;
  //set remaining sigma points
  for (int i = 0; i < n_x; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_ +n_x) * A.col(i);
    Xsig.col(i+1+n_x) = x_ - sqrt(lambda_ +n_x) * A.col(i);
  }

  //--------------------------------------------------------------------------------
  //Augmentation
  //--------------------------------------------------------------------------------

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;

  lambda_ = 3 - n_aug;
  for (int i = 0; i< n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug) * L.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda_+n_aug) * L.col(i);
  }
  //--------------------------------------------------------------------------------
  //Predict Sigma Points
  //--------------------------------------------------------------------------------


  //predict sigma points
  for (int i = 0; i< 2*n_aug+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) 
    {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else
    {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  //---------------------------------------------------------------------------------------
  //Predict mean and Covariance
  //---------------------------------------------------------------------------------------


  //create vector for weights
  lambda_ = 3 - n_aug;
 

 
  //set weights
  for (int i=0; i <(2*n_aug+1) ; i++)
  {
      if (i ==0)
      {
         weights_(i) = lambda_/(lambda_+n_aug);  
      }
      else
      {
        weights_(i) = 1/(2*(lambda_+n_aug)); 
      }
  }

  x_.fill(0);
  for (int j =0;j<(2*n_aug+1);j++)
  {
      x_ = x_ + weights_(j)*Xsig_pred_.col(j);
  }

  //predict state covariance matrix
    P_.fill(0.0);
    VectorXd x_diff = VectorXd(n_x);
    for (int j =0;j<(2*n_aug+1);j++)
    {
        x_diff = (Xsig_pred_.col(j)-x_);
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        P_ = P_ + weights_(j)*x_diff*x_diff.transpose();
    }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) 
{
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //---------------------------------------------------------------------------------------
  //Predict measurements
  //---------------------------------------------------------------------------------------
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;

//set vector for weights

  /*double weight_0 = lambda_/(lambda_+n_aug);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  
    double weight = 0.5/(n_aug+lambda_);
    weights_(i) = weight;
  }
  */

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_*std_laspx_,0,
        0,std_laspy_*std_laspy_;
  Zsig.fill(0.0);
  //transform sigma points into measurement space
  for (int j=0; j<(2 * n_aug + 1);j++)
  {
    Zsig(0,j) = Xsig_pred_(0,j); //px
    Zsig(1,j) = Xsig_pred_(1,j); //py
  }


  //calculate mean predicted measurement
  
    z_pred.fill(0.0);
    for (int i=0;i<(2*n_aug+1);i++)
    {
        z_pred = z_pred +  Zsig.col(i)*weights_(i);
    }
  //calculate innovation covariance matrix S

  S.fill(0.0);
  VectorXd z_diff = VectorXd(n_z);
  for (int i=0;i<(2 * n_aug + 1);i++)
  {
      z_diff = (Zsig.col(i) -z_pred );
      
      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
      S = S +  weights_(i)* z_diff *  z_diff.transpose();
  }
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  z_diff = meas_package.raw_measurements_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) 
{
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //---------------------------------------------------------------------------------------
  //Predict measurements
  //---------------------------------------------------------------------------------------
    //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda_/(lambda_+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) 
  {  
    double weight = 0.5/(n_aug+lambda_);
    weights(i) = weight;
  }

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);


  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_,0,0,
        0,std_radphi_*std_radphi_,0,
        0,0,std_radrd_*std_radrd_;

  //transform sigma points into measurement space
  for (int j=0; j<(2 * n_aug + 1);j++)
  {
    double px       =  Xsig_pred_(0,j);
    double py       =  Xsig_pred_(1,j);
    double v        =  Xsig_pred_(2,j);
    double phi      =  Xsig_pred_(3,j);
    //double phi_dot  =  Xsig_pred_(4,j);
    
    Zsig(0,j) = sqrt(px*px +py*py); //rho
    Zsig(1,j) = atan2(py,px); //phi
    Zsig(2,j) = (px*cos(phi)*v + py*sin(phi)*v)/(sqrt(px*px+py*py)); //rho_dot
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0;i<(2*n_aug+1);i++)
  {
      z_pred = z_pred +  Zsig.col(i)*weights(i);
  }
  
  //calculate innovation covariance matrix S

  S.fill(0.0);
  VectorXd z_diff = VectorXd(n_z);
  for (int i=0;i<(2 * n_aug + 1);i++)
  {
      z_diff = (Zsig.col(i) -z_pred );
      
      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
      S = S +  weights(i)* z_diff *  z_diff.transpose();
  }
  S = S + R;

    //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual

    z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  z_diff = meas_package.raw_measurements_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
