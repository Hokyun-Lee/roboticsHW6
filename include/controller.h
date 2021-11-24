#ifndef CONTROLLER_H_
#define CONTROLLER_H_

#include <iostream>
#include <Eigen/Dense>
#include <rbdl/rbdl.h>
#include <memory>
#include <fstream>
#include "math_type_define.h"
#include "franka_model_interface.h"

#define EYE(X) Matrix<double, X, X>::Identity()

using namespace RigidBodyDynamics;
using namespace std;
using namespace Eigen;
using namespace DyrosMath;

class ArmController
{
    size_t dof_;

	// Initial state
	Vector7d q_init_;
	Vector7d qdot_init_;

	// Current state
	Vector7d q_;
	Vector7d qdot_;
	Vector7d qddot_;
	Vector7d torque_;
	Vector7d q_error_sum_;

	// Control value (position controlled)
	Vector7d q_desired_; // Control value
	Vector7d torque_desired_;

	// Task space
	Vector3d x_init_;
	Vector3d x_init_2_;
	Vector3d x_dot_init_;
	Vector3d x_;
	Vector3d x_2;
	Matrix3d rotation_;
	Matrix3d rotation_init_;
	Matrix3d rotation_target_;
	Matrix3d rotation_cubic_;
	Vector3d rotation_cubic_dot;
	Vector3d phi_;
	Vector6d x_dot_; // 6D (linear + angular)
	Vector6d x_dot_2_; // 6D (linear + angular)
	Vector6d x_error_; 
    Vector3d x_from_q_desired_;
    Matrix3d rotation_from_q_desired_;

	// Dynamics
	Vector7d g_; // Gravity torque
	Matrix7d m_; // Mass matrix
	Matrix7d m_inverse_; // Inverse of mass matrix

	// For controller
	Matrix<double,3,7> j_v_;	// Linear velocity Jacobian matrix
	Matrix<double, 3, 7> j_w_;	// Angular veolicty Jacobain matrix
	Matrix<double, 6, 7> j_;	// Full basic Jacobian matrix
	Matrix<double, 7, 6> j_inverse_;	// Jacobain inverse storage 
	Matrix<double, 6, 7> j_2;
	Matrix<double, 3, 7> j_2_v_;
    Matrix<double, 6, 7> j_from_q_desired_;

	VectorXd q_temp_;	// For RBDL 
	VectorXd qdot_temp_;
	VectorXd qddot_temp_;
	MatrixXd j_temp_;	// For RBDL 
	MatrixXd j_temp_2;
	MatrixXd m_temp_;
	VectorXd g_temp_;   // For RBDL 

	Vector7d q_cubic_;
	Vector7d qdot_cubic_;
	Vector7d q_target_;
	Vector7d qdot_target_;

	Vector3d x_cubic_;
	Vector3d x_cubic_2_;
	Vector3d x_cubic_old_;
	Vector6d x_cubic_dot_;
	Vector3d x_cubic_dot_1_;
	Vector3d x_cubic_dot_2_;
	Vector3d x_target_;
	Vector3d x_target_2_;
	Vector3d x_dot_target_;

	double h1;
	double h2;

    unsigned long tick_;
    double play_time_;
    double hz_;
    double dt_;
    double control_start_time_;

    std::string control_mode_;
    bool is_mode_changed_;
	bool set_target_;

	Matrix<double, 6, 6> Kp_jacobian_;
	// HW5-1
	Vector7d target_q;
	Vector3d target_x;
	Vector3d target_x_save;
	Vector3d target_x_dot;
	Matrix3d target_ori;
	Vector3d x_dot;
	Vector3d ori_dot;
	Vector3d f_star;
	Vector3d m_star;
	Matrix3d Kp_, Kv_;
	Matrix7d Kp_joint_, Kv_joint_;
	Matrix6d Lamda;
	Vector6d control_input;
	Vector7d control_input_2;
	Matrix<double, 6, 7> j_t_dyn_cons_inv;

	//HW6-1
	Vector3d obs_x;
	Vector3d obs_dist;
	Vector3d f_star_obs;
	double k_obs;


	// file write
	string path;
	ofstream writeFile;
	bool write_;
	bool write_init_;

    FrankaModelInterface model_interface_;

private:
    void printState();
	//void moveJointPosition(const Vector7d &target_position, double duration);
	//void moveJointPositionTorque(const Vector7d &target_position, double duration);

public:
	void closeFile();
	void readData(const Vector7d &position, const Vector7d &velocity, const Vector7d &torque);
	void readData(const Vector7d &position, const Vector7d &velocity);
	const Vector7d & getDesiredPosition();
	const Vector7d & getDesiredTorque();

public:
		ArmController(double hz, franka::Model &model) :tick_(0), play_time_(0.0), hz_(hz),control_mode_("none"), is_mode_changed_(false),model_interface_(model)
	{
			initDimension();
	}


    void setMode(const std::string & mode);
    void initDimension();
    void initPosition(franka::RobotState state);
    void compute();
    void updateTime(double dt);
};

#endif
