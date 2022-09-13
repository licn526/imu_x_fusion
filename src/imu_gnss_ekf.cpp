#include <ros/ros.h>
#include <sensor_msgs/NavSatFix.h>

#include <deque>
#include <fstream>
#include <iostream>

#include "common/view.hpp"
#include "estimator/ekf.hpp"
#include "sensor/gnss.hpp"

Eigen::Vector3d sum_of_res(0, 0, 0), expect_of_res(0, 0, 0);
long long rescnt = 0;
double ekfla, ekflon, ekfu;
double zero_time = -1;

const double kDegreeToRadian = M_PI / 180.;
const double sigma_pv = 10;
const double sigma_rp = 10 * kDegreeToRadian;
const double sigma_yaw = sigma_rp * sigma_rp;
const double sigma_ba = 0.02;
const double sigma_bg = 0.02;

namespace cg {

ANGULAR_ERROR State::kAngError = ANGULAR_ERROR::LOCAL_ANGULAR_ERROR;

class FusionNode {
 public:
  FusionNode(ros::NodeHandle &nh) : viewer_(nh) {
    double acc_n, gyr_n, acc_w, gyr_w;
    nh.param("asdfsadf", acc_n, 1e-2);
    nh.param("adsfasdf", gyr_n, 1e-4);
    nh.param("asdfasdfasd", acc_w, 1e-6);
    nh.param("adsfasdfad", gyr_w, 1e-8);

    

    ekf_ptr_ = std::make_unique<EKF>(acc_n, gyr_n, acc_w, gyr_w);
    ekf_ptr_->state_ptr_->set_cov(sigma_pv, sigma_pv, sigma_rp, sigma_yaw, sigma_ba, sigma_bg);

    ekf_ptr_->observer_ptr_ = std::make_shared<GNSS>();

    std::string topic_imu = "/imu/data";
    std::string topic_gps = "/imu/nav_sat_fix";
    

    imu_sub_ = nh.subscribe<sensor_msgs::Imu>(topic_imu, 10, boost::bind(&EKF::imu_callback, ekf_ptr_.get(), _1));
    gps_sub_ = nh.subscribe<sensor_msgs::NavSatFix>(topic_gps, 10, &FusionNode::gps_callback, this);
    
    std::string topic_ekf = "/sbg/ekf_nav";
    ekf_sub_ = nh.subscribe<sensor_msgs::NavSatFix>(topic_ekf, 10, &FusionNode::ekf_callback, this);

    // log files
    observability.open("observability.csv");
    file_state_.open("fusion_state.csv");
    //file_imu_data.open("imu_data.csv")
  }

  ~FusionNode() {
    if (observability.is_open()) observability.close();
    if (file_state_.is_open()) file_state_.close();
    //if (file_imu_data.is_open()) file_imu_data.close();
  }
  

  void gps_callback(const sensor_msgs::NavSatFixConstPtr &gps_msg);
  void ekf_callback(const sensor_msgs::NavSatFixConstPtr &ekf_msg);
 private:
  ros::Subscriber imu_sub_;
  ros::Subscriber gps_sub_;
  ros::Subscriber ekf_sub_;

  EKFPtr ekf_ptr_;
  Viewer viewer_;

  std::ofstream observability;
  std::ofstream file_state_;
  //std::ofstream file_imu_data;
};

void FusionNode::ekf_callback(const sensor_msgs::NavSatFixConstPtr &ekf_msg){
    ekfla = ekf_msg->latitude;
    ekflon = ekf_msg->longitude;
    ekfu = ekf_msg->altitude;
    std::cout << "*\n";
}

void FusionNode::gps_callback(const sensor_msgs::NavSatFixConstPtr &gps_msg) {
  /*if (gps_msg->status.status != 2) {
    printf("[cggos %s] ERROR: Bad GPS Message!!!\n", __FUNCTION__);
    return;
  }*/
  GpsDataPtr gps_data_ptr = std::make_shared<GpsData>();
  gps_data_ptr->timestamp = gps_msg->header.stamp.toSec();
  if(zero_time < 0) zero_time = gps_msg->header.stamp.toSec();
  gps_data_ptr->lla[0] = gps_msg->latitude;
  gps_data_ptr->lla[1] = gps_msg->longitude;
  gps_data_ptr->lla[2] = gps_msg->altitude;
  gps_data_ptr->cov = Eigen::Map<const Eigen::Matrix3d>(gps_msg->position_covariance.data());
  /*gps_data_ptr->cov.setZero();
  gps_data_ptr->cov(0, 0) = 25;
  gps_data_ptr->cov(1, 1) = 25;
  gps_data_ptr->cov(2, 2) = 25;*/
  
  if (!ekf_ptr_->inited_) {
    if (!ekf_ptr_->init(gps_data_ptr->timestamp)) return;

    std::dynamic_pointer_cast<GNSS>(ekf_ptr_->observer_ptr_)->set_params(gps_data_ptr);

    printf("[cggos %s] System initialized.\n", __FUNCTION__);

    return;
  }

  std::cout << "---------------------" << std::endl;

  const Eigen::Isometry3d &Twb = ekf_ptr_->state_ptr_->pose();
  const auto &p_G_Gps = std::dynamic_pointer_cast<GNSS>(ekf_ptr_->observer_ptr_)->g2l(gps_data_ptr);

  const auto &residual = ekf_ptr_->observer_ptr_->measurement_residual(Twb.matrix(), p_G_Gps);

  std::cout << "number of res:" << ++rescnt << std::endl;
  std::cout << "latitude: " << gps_data_ptr->lla[0] << ' ' << "longitude: " << gps_data_ptr->lla[1] << ' ' << "altitude: " << gps_data_ptr->lla[2] << std::endl;
  std::cout << "res: " << residual.transpose() << std::endl;
  
  sum_of_res += residual;
  expect_of_res = sum_of_res / rescnt;

  std::cout << "sum of res:" << sum_of_res.transpose() << std::endl;
  std::cout << "expect of res: " << expect_of_res.transpose() << std::endl;


  const auto &H = ekf_ptr_->observer_ptr_->measurement_jacobian(Twb.matrix(), p_G_Gps);

  Eigen::Matrix<double, kStateDim, 3> K;

  const Eigen::Matrix3d &R = gps_data_ptr->cov;

  ekf_ptr_->update_K(H, R, K);
  ekf_ptr_->update_P(H, R, K);
  *ekf_ptr_->state_ptr_ = *ekf_ptr_->state_ptr_ + K * residual;

  std::cout << "acc bias: " << ekf_ptr_->state_ptr_->acc_bias.transpose() << std::endl;
  std::cout << "gyr bias: " << ekf_ptr_->state_ptr_->gyr_bias.transpose() << std::endl;
  // print P(15 * 15)
  std::cout << "P:\n";
  for (int i = 0; i < 15; ++i) std::cout << ekf_ptr_->state_ptr_->cov(i, i) << ' ';
  std::cout << std::endl; // 15 * 15
  // print K(15 * 3)
  std::cout << "K:\n";
  for (int i = 0; i < 15; ++i) std::cout << K(i, i % 3) << ' ';
  std::cout << std::endl;  // 15 * 3
  // print R(3 * 3)
  std::cout << "R:\n";
  std::cout << R(0, 0) << ' ' << R(1, 1) << ' ' << R(2, 2);
  std::cout << std::endl;  // 3 * 3
  // print x(15 * 1)
  std::cout << "x:\n"
            << ekf_ptr_->state_ptr_->p_wb_.transpose() << std::endl 
            << ekf_ptr_->state_ptr_->v_wb_.transpose() << std::endl
            << ekf_ptr_->state_ptr_->Rwb_.transpose() << std::endl
            << ekf_ptr_->state_ptr_->acc_bias.transpose() << std::endl
            << ekf_ptr_->state_ptr_->gyr_bias.transpose() << std::endl;
  // print H
  std::cout << H << std::endl;
  // std::cout << "Q:\n" << noise_cov_discret_time(dt) << std::endl;
  std::cout << "---------------------" << std::endl;

  // save data
  {
    viewer_.publish_gnss(*ekf_ptr_->state_ptr_);

    // save state p q lla
    const auto &lla = std::dynamic_pointer_cast<GNSS>(ekf_ptr_->observer_ptr_)->l2g(ekf_ptr_->state_ptr_->p_wb_);

    const Eigen::Quaterniond q_GI(ekf_ptr_->state_ptr_->Rwb_);
    file_state_ << std::fixed << std::setprecision(15) << ekf_ptr_->state_ptr_->timestamp << ", "
                 << ekf_ptr_->state_ptr_->p_wb_[0] << ", " << ekf_ptr_->state_ptr_->p_wb_[1] << ", "
                << ekf_ptr_->state_ptr_->p_wb_[2] << ", " << ekf_ptr_->state_ptr_->v_wb_[0] << ", "
                << ekf_ptr_->state_ptr_->v_wb_[1] << ", " << ekf_ptr_->state_ptr_->v_wb_[2] << ", " << q_GI.x() << ", "
                << q_GI.y() << ", " << q_GI.z() << ", "
                 << q_GI.w() << ", " << lla[0] << ", " << lla[1] << ", " << lla[2] << ", " << gps_data_ptr->lla[0]
                 << ", " << gps_data_ptr->lla[1] << ", " << gps_data_ptr->lla[2] << ", " << ekfla << ", " << ekflon
                 << ", " << ekfu
                 << std::endl;

    observability << std::fixed << std::setprecision(15) 
                  << gps_data_ptr->timestamp - zero_time << ", "
                  << sqrt(sigma_pv * sigma_pv / ekf_ptr_->state_ptr_->cov(0, 0)) << ", "
                  << sqrt(sigma_pv * sigma_pv / ekf_ptr_->state_ptr_->cov(1, 1)) << ", " 
                  << sqrt(sigma_pv * sigma_pv / ekf_ptr_->state_ptr_->cov(2, 2)) << ", " 
                  << sqrt(sigma_pv * sigma_pv / ekf_ptr_->state_ptr_->cov(3, 3)) << ", "
                  << sqrt(sigma_pv * sigma_pv / ekf_ptr_->state_ptr_->cov(4, 4)) << ", " 
                  << sqrt(sigma_pv * sigma_pv / ekf_ptr_->state_ptr_->cov(5, 5)) << ", " 
                  << sqrt(sigma_rp * sigma_rp / ekf_ptr_->state_ptr_->cov(6, 6)) << ", "
                  << sqrt(sigma_rp * sigma_rp / ekf_ptr_->state_ptr_->cov(7, 7)) << ", " 
                  << sqrt(sigma_yaw / ekf_ptr_->state_ptr_->cov(8, 8)) << ", "
                  << sqrt(sigma_ba * sigma_ba / ekf_ptr_->state_ptr_->cov(9, 9)) << ", "
                  << sqrt(sigma_ba * sigma_ba / ekf_ptr_->state_ptr_->cov(10, 10)) << ", "
                  << sqrt(sigma_ba * sigma_ba / ekf_ptr_->state_ptr_->cov(11, 11)) << ", "
                  << sqrt(sigma_bg * sigma_bg / ekf_ptr_->state_ptr_->cov(12, 12)) << ", "
                  << sqrt(sigma_bg * sigma_bg / ekf_ptr_->state_ptr_->cov(13, 13)) << ", "
                  << sqrt(sigma_bg * sigma_bg / ekf_ptr_->state_ptr_->cov(14, 14)) << std::endl;

     /*<< ekf_ptr_->state_ptr_->cov(0, 0) / (sigma_pv * sigma_pv) << ", " << ekf_ptr_->state_ptr_->cov(1, 1) / (sigma_pv * sigma_pv) << ", " << ekf_ptr_->state_ptr_->cov(2, 2) / (sigma_pv * sigma_pv) << ", " << ekf_ptr_->state_ptr_->cov(3, 3) / (sigma_pv * sigma_pv) << ", " << ekf_ptr_->state_ptr_->cov(4, 4) / (sigma_pv * sigma_pv) << ", " << ekf_ptr_->state_ptr_->cov(5, 5) / (sigma_pv * sigma_pv) << ", " << ekf_ptr_->state_ptr_->cov(6, 6) / (sigma_rp * sigma_rp) << ", " << ekf_ptr_->state_ptr_->cov(7, 7) / (sigma_rp * sigma_rp) << ", " << ekf_ptr_->state_ptr_->cov(8, 8) / sigma_yaw << ", " << ekf_ptr_->state_ptr_->cov(9, 9) / (sigma_ba * sigma_ba) << ", " << ekf_ptr_->state_ptr_->cov(10, 10) / (sigma_ba * sigma_ba) << ", " << ekf_ptr_->state_ptr_->cov(11, 11) / (sigma_ba * sigma_ba) << ", " << ekf_ptr_->state_ptr_->cov(12, 12) / (sigma_bg * sigma_bg) << ", " << ekf_ptr_->state_ptr_->cov(13, 13) / (sigma_bg * sigma_bg) << ", " << ekf_ptr_->state_ptr_->cov(14, 14) / (sigma_bg * sigma_bg) << std::endl;*/
    
  }
}

}  // namespace cg

int main(int argc, char **argv) {
  ros::init(argc, argv, "imu_gnss_fusion");

  ros::NodeHandle nh;
  cg::FusionNode fusion_node(nh);

  ros::spin();


  return 0;
}
