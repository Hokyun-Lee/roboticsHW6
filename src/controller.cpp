#include "controller.h"
#include <iostream>
#include <iomanip>
#include <cmath>


void ArmController::compute()
{
    // Kinematics and dynamics calculation ------------------------------
    Eigen::Affine3d transform = model_interface_.getTransform(franka::Frame::kEndEffector, q_);
    Eigen::Affine3d transform2 = model_interface_.getTransform(franka::Frame::kJoint4, q_);
    x_ = transform.translation();
    rotation_ = transform.linear();
    j_ = model_interface_.getJacobianMatrix(franka::Frame::kEndEffector, q_);
    m_ = model_interface_.getMassMatrix(q_);
    m_inverse_ = m_.inverse();
    // For CLIK
    Eigen::Affine3d transform_from_q_desired = model_interface_.getTransform(franka::Frame::kEndEffector,
    q_desired_);
    x_from_q_desired_ = transform_from_q_desired.translation();
    rotation_from_q_desired_ = transform_from_q_desired.linear();
    j_from_q_desired_ = model_interface_.getJacobianMatrix(franka::Frame::kEndEffector, q_desired_);
    // -----------------------------------------------------

    
	if (is_mode_changed_)
	{
		is_mode_changed_ = false;
		write_ = false;
		write_init_ = true;

		control_start_time_ = play_time_;

		q_init_ = q_;
		qdot_init_ = qdot_;
		q_error_sum_.setZero();
		q_desired_ = q_;

		x_init_ = x_;
		x_init_2_ = x_2;
		x_dot_init_ = x_dot_.block<3,1>(0,0);
		x_cubic_old_ = x_;
		rotation_init_ = rotation_;
	}
	

	if (control_mode_ == "joint_ctrl_home")
	{
        double duration = 5.0;
		Vector7d target_position;
		target_position << 0.0, 0.0, 0.0, -M_PI / 4, 0.0, M_PI / 2, 0;
		q_desired_ = cubicVector<7>(play_time_, control_start_time_, control_start_time_ + duration, q_init_, target_position, qdot_init_, qdot_target_);
	}
	else if(control_mode_ == "joint_ctrl_init")
	{
        double duration = 5.0;
		Vector7d target_position;
		target_position << 0.0, 0.0, 0.0, -M_PI / 2., 0.0, M_PI / 2, 0;
		q_desired_ = cubicVector<7>(play_time_, control_start_time_, control_start_time_ + duration, q_init_, target_position, qdot_init_, qdot_target_);
	}
	else if (control_mode_ == "CLIK")
	{
		x_target_ << 0.25, 0.28, 0.65;
		rotation_target_ << 0.7071, 0.7071, 0, 0.7071, -0.7071, 0, 0, 0,-1;
        
		double duration = 3.0;

		Matrix<double, 7, 6> pseudo_j_;
		pseudo_j_ = j_from_q_desired_.transpose() * (j_from_q_desired_ * j_from_q_desired_.transpose()).inverse();
		// cout << "pseudo_j_ :\n" << pseudo_j_ << endl;

		for (int i = 0; i < 3; i++)
		{
			x_cubic_(i) = DyrosMath::cubicDot(play_time_, control_start_time_, control_start_time_ + duration, x_init_(i), x_target_(i), 0.0, 0.0, hz_);
		}

		xdot_desired_(0) = x_cubic_(0);
		xdot_desired_(1) = x_cubic_(1);
		xdot_desired_(2) = x_cubic_(2);
		xdot_desired_(3) = -DyrosMath::getPhi(rotation_, body_to_ee_rotation)(0);
		xdot_desired_(4) = -DyrosMath::getPhi(rotation_, body_to_ee_rotation)(1);
		xdot_desired_(5) = -DyrosMath::getPhi(rotation_, body_to_ee_rotation)(2);
		x_o_init_(0) = x_init_(0);
		x_o_init_(1) = x_init_(1);
		x_o_init_(2) = x_init_(2);
		x_o_init_(3) = rotation_init_(0);
		x_o_init_(4) = rotation_init_(1);
		x_o_init_(5) = rotation_init_(2);
		// cout << "xdot_desired_ :\n" << xdot_desired_ << endl;

		if (control_start_time_ == play_time_) {
			x_desired_ = x_o_init_ + xdot_desired_ * dt_;
			// cout << "ONCE" << endl;
		}
		else {
			x_desired_ = x_desired_ + xdot_desired_ * dt_;
		}

		Kp.setIdentity();
		Kp = 1 * Kp;

		x_q_desired_(0) = x_from_q_desired_(0);
		x_q_desired_(1) = x_from_q_desired_(1);
		x_q_desired_(2) = x_from_q_desired_(2);
		x_q_desired_(3) = rotation_from_q_desired_(0);
		x_q_desired_(4) = rotation_from_q_desired_(1);
		x_q_desired_(5) = rotation_from_q_desired_(2);

		//cout << "xdesired error\n" << x_desired_ - x_q_desired_ << endl;
		qdot_desired_ = pseudo_j_ * (xdot_desired_ + Kp * (x_desired_ - x_q_desired_));
		q_desired_ = q_desired_ + qdot_desired_ * 1 * dt_;
	}
    else if (control_mode_ == "Dynamics")
    {
        q_target_ << 0.0, 0.0, 0.0, -60*DEG2RAD, 0.0, 90*DEG2RAD, 0.0;
		
        // torque_desired_ = 
    }
	else if (control_mode_ == "VelocitySaturation")
	{
		double duration = 5.0;

		Lamda = (j_ * m_.inverse() * j_.transpose()).inverse();
		Matrix<double, 7, 6> j_bar_;
		j_bar_ = m_.inverse() * j_.transpose() * Lamda;

		target_x(0) = 0.3;
		target_x(1) = -0.012;
		target_x(2) = 0.52;

		Vector3d x_desired_dot_;
		x_desired_dot_ << 0, 0, 0;

		for (int i = 0; i < 3; i++) {
			x_desired_(i) = target_x(i);
			x_desired_dot_(i) = 0;
		}

		double Kp_val, Kv_val;

		Kp_val = 400.0;
		Kv_val = 40.0;

		Vector3d saturation_x_desired_dot;
		for (int i = 0; i < 3; i++) {
			saturation_x_desired_dot(i) = Kp_val / Kv_val * (x_desired_(i) - x_(i));
			if (abs(saturation_x_desired_dot(i)) > 0.3) {
				saturation_x_desired_dot(i) = 0.3 / abs(x_desired_(i) - x_(i)) * (x_desired_(i) - x_(i));
			}
		}

		f_star = Kv_val * EYE(3) * (saturation_x_desired_dot - x_dot_.head(3));

		m_star = -Kp_val * EYE(3) * DyrosMath::getPhi(rotation_, rotation_init_) - Kv_val * x_dot_.tail<3>();

		Vector6d F_zero_star_;
		F_zero_star_ << f_star, m_star;

		Vector7d torque_zero_;
		torque_zero_ = m_ * (Kp_val * EYE(7) * (q_init_ - q_) - Kv_val * EYE(7) * qdot_);

		torque_desired_ = j_.transpose() * Lamda * F_zero_star_ + (EYE(7) - j_.transpose() * j_bar_.transpose()) * torque_zero_ + g_;
	}
	else 
	{
		torque_desired_.setZero();
	}

	printState();

	tick_++;
}


void ArmController::printState()
{
	// TODO: Modify this method to debug your code

	static int DBG_CNT = 0;
	if (DBG_CNT++ > hz_ / 20.)
	{
		DBG_CNT = 0;

		//cout << "q desired:\t";
		//cout << std::fixed << std::setprecision(3) << q_desired_.transpose() << endl;
		//cout << "q cubic:\t";
		//cout << std::fixed << std::setprecision(3) << q_cubic_.transpose() << endl;
		//cout << "q now    :\t";
		//cout << std::fixed << std::setprecision(3) << q_.transpose() << endl;
		//cout << "x desired:\t";
		//cout << x_cubic_.transpose() << endl;
		//cout << "x        :\t";
		//cout << x_.transpose() << endl;

	}

	if (play_time_ == control_start_time_ + 3)
	{
		cout << "---------------------------------------------------------------------" << endl;
		cout << "                     control time finished                           " << endl;
		cout << "---------------------------------------------------------------------" << endl;
	}
}



// Controller Core Methods ----------------------------

void ArmController::setMode(const std::string & mode)
{
	is_mode_changed_ = true;
	control_mode_ = mode;
	cout << "Current mode (changed) : " << mode << endl;
}
void ArmController::initDimension()
{
	dof_ = DOF;
	q_temp_.resize(DOF);
	j_temp_.resize(6, DOF);
	j_temp_2.resize(6, DOF);
	j_temp_2.setZero();

	qddot_.setZero();

	x_target_.setZero();
	q_desired_.setZero();
	qdot_target_.setZero();
	torque_desired_.setZero();
	x_dot_init_.setZero();
	x_dot_target_.setZero();

	g_temp_.resize(DOF);
	m_temp_.resize(DOF, DOF);

	Kp_.setZero();
	Kv_.setZero();
	Kp_joint_.setZero();
	Kv_joint_.setZero();
}

void ArmController::readData(const Vector7d &position, const Vector7d &velocity, const Vector7d &torque)
{
	for (size_t i = 0; i < dof_; i++)
	{
		q_(i) = position(i);
		qdot_(i) = velocity(i);
		torque_(i) = torque(i);
	}
}
void ArmController::readData(const Vector7d &position, const Vector7d &velocity)
{
	for (size_t i = 0; i < dof_; i++)
	{
		q_(i) = position(i);
		qdot_(i) = velocity(i);
		torque_(i) = 0;
	}
}

const Vector7d & ArmController::getDesiredPosition()
{
	return q_desired_;
}

const Vector7d & ArmController::getDesiredTorque()
{
	return torque_desired_;
}

void ArmController::closeFile()
{
	if(writeFile.is_open())
		writeFile.close();
}

void ArmController::initPosition(franka::RobotState state)
{
    q_init_ = q_;
    q_desired_ = q_init_;
    model_interface_.setRobotState(state);
}

void ArmController::updateTime(double dt)
{
    dt_ = dt;
    play_time_ += dt;
}
// ----------------------------------------------------

