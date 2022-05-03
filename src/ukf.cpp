#include "ukf.h"
#include "Eigen/Dense"
#include "iostream"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  // Don't change sensor noise values:
  // ##################################################################
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
  // ##################################################################



  n_x_ = 5;
  n_aug_ = 7;
  Xsig_pred_ = Eigen::MatrixXd(n_x_, 2 * n_aug_ + 1);

  lambda_ = 3 - n_aug_;

  weights_ = Eigen::VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // std::cout << "Weights: " << weights_ << std::endl;

  bool is_weight_dim_ok = weights_.size() == 2*n_aug_+1;
  assert(is_weight_dim_ok && "Weights vector dimension is not ok!");


  // Define measurement noise covariance matrices for radar and laser:
  // Stochastic noise for radar measurements
  Noise_radar_ = Eigen::MatrixXd::Zero(3, 3);
  Noise_radar_ <<
    square(std_radr_), 0, 0,
    0, square(std_radphi_), 0,
    0, 0, square(std_radrd_);
  // Stochastic noise for laser measurements
  Noise_laser_ = Eigen::MatrixXd::Zero(2, 2);
  Noise_laser_ <<
    square(std_laspx_), 0,
    0, square(std_laspy_);

  use_laser_ = true;
  use_radar_ = true;

  // It will be updated first iteration of ProcessMeasurement():
  is_initialized_ = false;
  // Tuning:
  std_a_ = 3.5;
  std_yawdd_ = 1;
}

UKF::~UKF() {}

void
UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  // Initialization first iter:
  if (!is_initialized_)
  {
    // use both radar and laser or only radar initialization:
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_d = meas_package.raw_measurements_(2);
      // Use formulas from course:
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      // Initialize velocities and magnitude of velocity:
      double vx = rho_d * cos(phi);
      double vy = rho_d * sin(phi);
      double v = sqrt(vx * vx + vy * vy); // assumption: vx and vy are not zero initially
      // Initialize state vector:
      x_ << px, py, v, 0, 0;
      // Initialize covariance matrix:
      P_ <<
        square(std_radr_), 0, 0, 0, 0,
        0, square(std_radr_), 0, 0, 0,
        0, 0, square(std_radrd_), 0, 0,
        0, 0, 0, square(std_radphi_), 0,
        0, 0, 0, 0, square(std_radphi_);
    }
    // only laser initialization:
    // except the position, others are hidden
    else
    {
      x_ <<
          meas_package.raw_measurements_(0),
          meas_package.raw_measurements_(1),
          0., 0, 0;
      P_ <<
          square(std_laspx_), 0, 0, 0, 0,
          0, square(std_laspy_), 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;

    std::cout << "Initialization completed!" << std::endl;
    return;
  }
  // Time difference between current and previous measurements:
  double dt = (meas_package.timestamp_ - time_us_)/1e6; // in seconds
  time_us_ = meas_package.timestamp_;
  // First predict, then update:
  Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    UpdateRadar(meas_package);
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    UpdateLidar(meas_package);
}




void
UKF::Prediction(double delta_t)
{
  Eigen::VectorXd x_aug = Eigen::VectorXd(n_aug_);
  Eigen::MatrixXd P_aug = Eigen::MatrixXd::Zero(n_aug_,n_aug_);
  Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd::Zero(n_aug_,2*n_aug_+1);

  // augmented cov matrix
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5, 5) = square(std_a_);
  P_aug(6, 6) = square(std_yawdd_);

  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  Eigen::MatrixXd mat_tmp = P_aug.llt().matrixL();
  //std::cout << mat_tmp << std::endl;
  //std::cout << "mat_tmp.rows() = " << mat_tmp.rows() << std::endl;
  //std::cout << "mat_tmp.cols() = " << mat_tmp.cols() << std::endl;
  //std::cout << "-----------------------------------------" << std::endl;

  // augmented sigma points:
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i<n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug+sqrt(lambda_+n_aug_)*mat_tmp.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug-sqrt(lambda_+ n_aug_)*mat_tmp.col(i);
  }


  // predict sigma points
  for (int i = 0; i < 2*n_aug_ +1; i++)
  {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_p;
    double py_p;
    if (fabs(yawd) > 0.001) // dont take nan values :)
    {
      // Formulas for solving integral :)
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v; // because constant velocity model
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd; // because constant yaw rate model

    // Noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // Update sigma points
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
  // Predict state vector
  x_.fill(0.0);
  for (int i = 0; i<2*n_aug_+1; i++)
    x_ = x_ + weights_(i)*Xsig_pred_.col(i);
  // Predict covariance matrix
  P_.fill(0.0);
  for (int i = 0; i <2*n_aug_+1; i++)
  {
    Eigen::VectorXd x_diff = Xsig_pred_.col(i)-x_;
    normalize_angle(x_diff,3);
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}






void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  int meas_dim_laser = 2;

  Eigen::VectorXd laser_meas = Eigen::VectorXd(meas_dim_laser);
  laser_meas << meas_package.raw_measurements_(0),
      meas_package.raw_measurements_(1);
  Eigen::MatrixXd laser_meas_sig = Eigen::MatrixXd(meas_dim_laser,
                                            2 * n_aug_ + 1);

  // mean laser predicted measurement
  Eigen::VectorXd laser_pred = Eigen::VectorXd(meas_dim_laser);
  laser_pred.fill(0.0);

  // meas cov laser:
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(
      meas_dim_laser, meas_dim_laser);


  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    // measurement model for laser, it measures px and py:
    laser_meas_sig(0, i) = Xsig_pred_(0, i);
    laser_meas_sig(1, i) = Xsig_pred_(1, i);
    // mean predicted measurement
    laser_pred += weights_(i) * laser_meas_sig.col(i);
  }

  for (int i = 0; i<2*n_aug_+1; i++)
  {
    Eigen::VectorXd res = laser_meas_sig.col(i)-laser_pred;
    // update meas cov laser:
    S += weights_(i) * res * res.transpose();
  }
  S += Noise_laser_;

  Eigen::MatrixXd cross_corr = Eigen::MatrixXd::Zero(
      n_x_, meas_dim_laser);

  for (int i = 0; i<2*n_aug_+1; i++)
  {
    // residual
    Eigen::VectorXd res_meas = laser_meas_sig.col(i) - laser_pred;
    // state difference
    Eigen::VectorXd diff = Xsig_pred_.col(i) - x_;
    // normalize angles
    normalize_angle(diff,3);
    cross_corr += weights_(i) * diff * res_meas.transpose();
  }
  // Kalman gain K;
  Eigen::MatrixXd K = cross_corr * S.inverse();
  // Residual
  Eigen::VectorXd meas_diff = laser_meas - laser_pred;

  // update state mean and covariance matrix
  x_ += K * meas_diff;
  P_ -= K * S * K.transpose();
  // NIS Laser:
  NIS_laser_ = meas_diff.transpose() * S.inverse() * meas_diff;
  std::cout << "NIS_laser: " << NIS_laser_ << std::endl;
}










void
UKF::UpdateRadar(MeasurementPackage meas_package)
{
  int meas_dim_radar = 3;
  Eigen::VectorXd radar_meas = Eigen::VectorXd(meas_dim_radar);
  double meas_rho = meas_package.raw_measurements_(0);
  double meas_phi = meas_package.raw_measurements_(1);
  double meas_rhod = meas_package.raw_measurements_(2);
  radar_meas << meas_rho, meas_phi, meas_rhod;

  Eigen::MatrixXd radar_meas_sig = Eigen::MatrixXd(
      meas_dim_radar, 2 * n_aug_ + 1);

  Eigen::VectorXd radar_pred = Eigen::VectorXd(
      meas_dim_radar);
  radar_pred.fill(0.0);

  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(meas_dim_radar, meas_dim_radar);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    radar_meas_sig(0, i) = sqrt(square(p_x) + square(p_y));
    radar_meas_sig(1, i) = atan2(p_y, p_x);
    radar_meas_sig(2, i) = (p_x * v1 + p_y * v2) / sqrt(square(p_x) + square(p_y));

    //calculate mean predicted measurement
    radar_pred += weights_(i) * radar_meas_sig.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    Eigen::VectorXd diff = radar_meas_sig.col(i) - radar_pred;
    normalize_angle(diff, 1);
    S += weights_(i) * diff * diff.transpose();
  }
  S += Noise_radar_;

  Eigen::MatrixXd cross_corr = Eigen::MatrixXd::Zero(n_x_, meas_dim_radar);

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    Eigen::VectorXd meas_diff = radar_meas_sig.col(i) - radar_pred;
    normalize_angle(meas_diff, 1);
    // state difference
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    normalize_angle(x_diff, 3);
    cross_corr += weights_(i) * x_diff * meas_diff.transpose();
  }

  //  Kalman gain K;
  Eigen::MatrixXd K = cross_corr * S.inverse();
  Eigen::VectorXd meas_diff = radar_meas - radar_pred;

  // angle normalization
  normalize_angle(meas_diff, 1);

  // update state mean and covariance matrix
  x_ += K * meas_diff;
  P_ -= K * S * K.transpose();

  // NIS radar
  NIS_radar_ = meas_diff.transpose() * S.inverse() * meas_diff;
  std::cout << "NIS_radar: " << NIS_radar_ << std::endl;
}


double
UKF::square(const double &x) {return x*x;}

void
UKF::normalize_angle(VectorXd &vec, const int &index)
{
  while (vec(index) > M_PI) vec(index) -= 2. * M_PI;
  while (vec(index) < -M_PI) vec(index) += 2. * M_PI;
}
