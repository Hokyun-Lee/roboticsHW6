#include "controller.h"
#include <iostream>
#include <iomanip>
#include <cmath>

FILE* fp1 = fopen("Jacobian.txt", "w");
FILE* fp2 = fopen("CLIK.txt", "w");
FILE* fp3 = fopen("Weighted.txt", "w");


void ArmController::compute()
{
	// ---------------------------------
	//
	// q_		: joint position
	// qdot_	: joint velocity
	// x_		: end-effector position 
	// j_		: end-effector basic jacobian
	// m_		: mass matrix
	//
	//-------------------------------------------------------------------

	// Kinematics and dynamics calculation
	q_temp_ = q_;
	qdot_temp_ = qdot_;

	RigidBodyDynamics::UpdateKinematicsCustom(*model_, &q_temp_, &qdot_temp_, NULL);
	x_ = CalcBodyToBaseCoordinates(*model_, q_, body_id_[DOF - 1], com_position_[DOF - 1], true);
	rotation_ = CalcBodyWorldOrientation(*model_, q_, body_id_[DOF - 1], true).transpose();
	Matrix3d body_to_ee_rotation;
	body_to_ee_rotation.setIdentity();
	body_to_ee_rotation(1, 1) = -1;
	body_to_ee_rotation(2, 2) = -1;
	rotation_ = rotation_ * body_to_ee_rotation; // To Match RBDL model and CoppeliaSim model
	CalcPointJacobian6D(*model_, q_, body_id_[DOF - 1], com_position_[DOF - 1], j_temp_, true);

	for (int i = 0; i<2; i++)
		j_.block<3, DOF>(i * 3, 0) = j_temp_.block<3, DOF>(3 - i * 3, 0);

	j_v_ = j_.block < 3, DOF>(0, 0);
	x_dot_ = j_ * qdot_;

	NonlinearEffects(*model_, q_, Vector7d::Zero(), g_temp_);
	CompositeRigidBodyAlgorithm(*model_, q_, m_temp_, true);

	g_ = g_temp_;
	m_ = m_temp_;
	m_inverse_ = m_.inverse();

	// Kinematics calculation using q_desired_ (CLIK)
	MatrixXd j_temp_from_q_desired_;
	j_temp_from_q_desired_.resize(6, DOF);
	j_temp_from_q_desired_.setZero();

	x_from_q_desired_ = CalcBodyToBaseCoordinates(*model_, q_desired_, body_id_[DOF - 1], com_position_[DOF - 1], false);

	rotation_from_q_desired_ = CalcBodyWorldOrientation(*model_, q_desired_, body_id_[DOF - 1], false).transpose();
	rotation_from_q_desired_ = rotation_from_q_desired_ * body_to_ee_rotation; // To Match RBDL model and CoppeliaSim model

	CalcPointJacobian6D(*model_, q_desired_, body_id_[DOF - 1], com_position_[DOF - 1], j_temp_from_q_desired_, false);

	for (int i = 0; i<2; i++)
		j_from_q_desired_.block<3, DOF>(i * 3, 0) = j_temp_from_q_desired_.block<3, DOF>(3 - i * 3, 0);


	// Control
	if (is_mode_changed_)
	{
		is_mode_changed_ = false;

		control_start_time_ = play_time_;

		q_init_ = q_;
		qdot_init_ = qdot_;
		q_error_sum_.setZero();

		x_init_ = x_;
		x_cubic_old_ = x_;
		rotation_init_ = rotation_;
	}

	if (control_mode_ == "joint_ctrl_home")
	{
		Vector7d target_position;
		target_position << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, M_PI / 4;
		double duration = 3.0;

		for (int i = 0; i < dof_; i++)
		{
			q_desired_(i) = DyrosMath::cubic(play_time_, control_start_time_, control_start_time_ + duration, q_init_(i), target_position(i), 0.0, 0.0);
		}
	}
	else if(control_mode_ == "joint_ctrl_init")
	{
		Vector7d target_position;
		target_position << 0.0, 0.0, 0.0, -M_PI / 2., 0.0, M_PI / 2, 0.0 ;
		double duration = 1.0;

		for (int i = 0; i < dof_; i++)
		{
			q_desired_(i) = DyrosMath::cubic(play_time_, control_start_time_, control_start_time_ + duration, q_init_(i), target_position(i), 0.0, 0.0);
		}
	}
	else if (control_mode_ == "torque_ctrl_dynamic")
	{
		// torque_desired_ = 
	}
	else if (control_mode_ == "hw2-1")
	{
		double duration = 3.0;

		Matrix<double, 7, 6> pseudo_j_;
		pseudo_j_ = j_from_q_desired_.transpose() * (j_from_q_desired_ * j_from_q_desired_.transpose()).inverse();
		// cout << "pseudo_j_ :\n" << pseudo_j_ << endl;

		x_target_ << 0.25, 0.28, 0.65;
		for (int i = 0; i < 3; i++)
		{
			x_cubic_(i) = DyrosMath::cubicDot(play_time_, control_start_time_, control_start_time_ + duration, x_init_(i), x_target_(i), 0.0, 0.0, hz_);
		}

		xdot_desired_(0) = x_cubic_(0);
		xdot_desired_(1) = x_cubic_(1);
		xdot_desired_(2) = x_cubic_(2);
		xdot_desired_(3) = 0;
		xdot_desired_(4) = 0;
		xdot_desired_(5) = 0;
		// cout << "xdot_desired_ :\n" << xdot_desired_ << endl;
		
		qdot_desired_ = pseudo_j_ * xdot_desired_;
		q_desired_ = q_desired_ + qdot_desired_ * 1 / hz_;

		fprintf(fp1, "%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t \n", play_time_-control_start_time_, x_(0), x_(1), x_(2), x_target_(0), x_target_(1), x_target_(2), q_desired_(0), q_desired_(1), q_desired_(2), q_desired_(3), q_desired_(4), q_desired_(5), q_desired_(6), q_(0), q_(1), q_(2), q_(3), q_(4), q_(5), q_(6), qdot_(3));
	}
	else if (control_mode_ == "hw2-2")
	{
		double duration = 3.0;

		Matrix<double, 7, 6> pseudo_j_;
		pseudo_j_ = j_from_q_desired_.transpose() * (j_from_q_desired_ * j_from_q_desired_.transpose()).inverse();
		// cout << "pseudo_j_ :\n" << pseudo_j_ << endl;

		x_target_ << 0.25, 0.28, 0.65;
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
			x_desired_ = x_o_init_ + xdot_desired_ * 1 / hz_;
			// cout << "ONCE" << endl;
		}
		else {
			x_desired_ = x_desired_ + xdot_desired_ * 1 / hz_;
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
		q_desired_ = q_desired_ + qdot_desired_ * 1 / hz_;

		fprintf(fp2, "%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t \n", play_time_ - control_start_time_, x_(0), x_(1), x_(2), x_target_(0), x_target_(1), x_target_(2), q_desired_(0), q_desired_(1), q_desired_(2), q_desired_(3), q_desired_(4), q_desired_(5), q_desired_(6), q_(0), q_(1), q_(2), q_(3), q_(4), q_(5), q_(6), qdot_(3));
	}
	else if (control_mode_ == "hw2-3")
	{
	    double duration = 3.0;

		// cout << "pseudo_j_ :\n" << pseudo_j_ << endl;

		x_target_ << 0.25, 0.28, 0.65;
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
			x_desired_ = x_o_init_ + xdot_desired_ * 1 / hz_;
			// cout << "ONCE" << endl;
		}
		else {
			x_desired_ = x_desired_ + xdot_desired_ * 1 / hz_;
		}
		
		Kp.setIdentity();
		Kp = 1 * Kp;

		x_q_desired_(0) = x_from_q_desired_(0);
		x_q_desired_(1) = x_from_q_desired_(1);
		x_q_desired_(2) = x_from_q_desired_(2);
		x_q_desired_(3) = rotation_from_q_desired_(0);
		x_q_desired_(4) = rotation_from_q_desired_(1);
		x_q_desired_(5) = rotation_from_q_desired_(2);

		W_inv_.setIdentity();
		W_inv_(3, 3) = 0.01;
		// W_inv_(3, 3) = 0.001;

		Matrix<double, 7, 6> pseudo_j_;
		pseudo_j_ = W_inv_ * j_from_q_desired_.transpose() * (j_from_q_desired_ * W_inv_ * j_from_q_desired_.transpose()).inverse();

		//cout << "xdesired error\n" << x_desired_ - x_q_desired_ << endl;
		qdot_desired_ = pseudo_j_ * (xdot_desired_ + Kp * (x_desired_ - x_q_desired_));
		q_desired_ = q_desired_ + qdot_desired_ * 1 / hz_;

		fprintf(fp3, "%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t \n", play_time_ - control_start_time_, x_(0), x_(1), x_(2), x_target_(0), x_target_(1), x_target_(2), q_desired_(0), q_desired_(1), q_desired_(2), q_desired_(3), q_desired_(4), q_desired_(5), q_desired_(6), q_(0), q_(1), q_(2), q_(3), q_(4), q_(5), q_(6), qdot_(3));
	}
	else
	{
		torque_desired_ = g_;
	}

	printState();

	tick_++;
	play_time_ = tick_ / hz_;	// second
}


void ArmController::printState()
{
	// TODO: Modify this method to debug your code

	static int DBG_CNT = 0;
	if (DBG_CNT++ > hz_ / 20.)
	{
		DBG_CNT = 0;

		cout << "q now    :\t";
		cout << std::fixed << std::setprecision(3) << q_.transpose() << endl;
		cout << "q desired:\t";
		cout << std::fixed << std::setprecision(3) << q_desired_.transpose() << endl;
		cout << "torque   :\t";
		cout << std::fixed << std::setprecision(3) << torque_.transpose() << endl;
		cout << "torque desired:\t";
		cout << std::fixed << std::setprecision(3) << torque_desired_.transpose() << endl;
		cout << "x        :\t";
		cout << x_.transpose() << endl;
		cout << endl;
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

	qddot_.setZero();

	x_target_.setZero();
	q_desired_.setZero();
	torque_desired_.setZero();

	g_temp_.resize(DOF);
	m_temp_.resize(DOF, DOF);
}

void ArmController::initModel()
{
    model_ = make_shared<Model>();

    model_->gravity = Vector3d(0., 0, -GRAVITY);

    double mass[DOF];
    mass[0] = 1.0;
    mass[1] = 1.0;
    mass[2] = 1.0;
    mass[3] = 1.0;
    mass[4] = 1.0;
    mass[5] = 1.0;
    mass[6] = 1.0;

    Vector3d axis[DOF];
	axis[0] = Eigen::Vector3d::UnitZ();
	axis[1] = Eigen::Vector3d::UnitY();
	axis[2] = Eigen::Vector3d::UnitZ();
	axis[3] = -1.0*Eigen::Vector3d::UnitY();
	axis[4] = Eigen::Vector3d::UnitZ();
	axis[5] = -1.0*Eigen::Vector3d::UnitY();
	axis[6] = -1.0*Eigen::Vector3d::UnitZ();


	Eigen::Vector3d global_joint_position[DOF];

	global_joint_position[0] = Eigen::Vector3d(0.0, 0.0, 0.3330);
	global_joint_position[1] = global_joint_position[0];
	global_joint_position[2] = Eigen::Vector3d(0.0, 0.0, 0.6490);
	global_joint_position[3] = Eigen::Vector3d(0.0825, 0.0, 0.6490);
	global_joint_position[4] = Eigen::Vector3d(0.0, 0.0, 1.0330);
	global_joint_position[5] = Eigen::Vector3d(0.0, 0.0, 1.0330);
	global_joint_position[6] = Eigen::Vector3d(0.0880, 0.0, 1.0330);

	joint_posision_[0] = global_joint_position[0];
	for (int i = 1; i < DOF; i++)
		joint_posision_[i] = global_joint_position[i] - global_joint_position[i - 1];

	com_position_[0] = Vector3d(0.000096, -0.0346, 0.2575);
	com_position_[1] = Vector3d(0.0002, 0.0344, 0.4094);
	com_position_[2] = Vector3d(0.0334, 0.0266, 0.6076);
	com_position_[3] = Vector3d(0.0331, -0.0266, 0.6914);
	com_position_[4] = Vector3d(0.0013, 0.0423, 0.9243);
	com_position_[5] = Vector3d(0.0421, -0.0103, 1.0482);
	com_position_[6] = Vector3d(0.1, -0.0120, 0.9536);

	for (int i = 0; i < DOF; i++)
		com_position_[i] -= global_joint_position[i];

    Math::Vector3d inertia[DOF];
	for (int i = 0; i < DOF; i++)
		inertia[i] = Eigen::Vector3d::Identity() * 0.001;

    for (int i = 0; i < DOF; i++) {
        body_[i] = Body(mass[i], com_position_[i], inertia[i]);
        joint_[i] = Joint(JointTypeRevolute, axis[i]);
        if (i == 0)
            body_id_[i] = model_->AddBody(0, Math::Xtrans(joint_posision_[i]), joint_[i], body_[i]);
        else
            body_id_[i] = model_->AddBody(body_id_[i - 1], Math::Xtrans(joint_posision_[i]), joint_[i], body_[i]);
    }
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



void ArmController::initPosition()
{
    q_init_ = q_;
    q_desired_ = q_init_;
}

// ----------------------------------------------------

