#pragma once
#include <iostream>
#include <string>
#include <functional>
#include <Eigen/Dense>


using namespace std;

extern "C" {
#include "extApi.h"
}

const std::string JOINT_HANDLE_PREFIX{ "panda_joint" };

class CoppeliaBridge
{
public:
	enum ControlMode { CTRL_POSITION, CTRL_VELOCITY, CTRL_TORQUE };

	CoppeliaBridge(ControlMode mode = CTRL_POSITION);
	~CoppeliaBridge();
	
	bool simConnectionCheck();
	void simLoop();

	void write();
	void read();

	void setDesiredPosition(const Eigen::Matrix<double, DOF, 1> & desired_q);
	void setDesiredTorque(const Eigen::Matrix<double, DOF, 1> & desired_torque);
	const Eigen::Matrix<double, DOF, 1> & getPosition();
	const Eigen::Matrix<double, DOF, 1> & getVelocity();

	const size_t getTick() { return tick_; }

private:
	Eigen::Matrix<double, DOF, 1> current_q_;
	Eigen::Matrix<double, DOF, 1> current_q_dot_;
	Eigen::Matrix<double, DOF, 1> desired_q_;
	Eigen::Matrix<double, DOF, 1> desired_torque_;

	simxInt clientID_;
	simxInt motorHandle_[DOF];	/// < Depends on simulation envrionment
	simxInt objectHandle_;

	size_t tick_{ 0 };
	
	ControlMode control_mode_;

	void simxErrorCheck(simxInt error);
	void simInit();
	void getHandle(); 	/// < Depends on simulation envrionment
};
