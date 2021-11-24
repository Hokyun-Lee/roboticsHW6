#include "controller.h"
#include <iostream>
#include <iomanip>
#include <cmath>

//FILE* fp1 = fopen("4-1.txt", "w");
//FILE* fp2 = fopen("4-2-1.txt", "w");
//FILE* fp3 = fopen("4-2-2.txt", "w");
//FILE* fp4 = fopen("4-3-1.txt", "w");
//FILE* fp5 = fopen("4-3-2.txt", "w");


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
	x2_ = CalcBodyToBaseCoordinates(*model_, q_, body_id_[DOF - 4], com_position_[DOF - 4], false);
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

	Kp = 400 * EYE(7);
	Kv = 40 * EYE(7);

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
		Vector7d target_q_;
		target_q_ << 0.0, 0.0, 0.0, - M_PI / 2, 0.0, M_PI / 2, 0.0;
		double duration = 1.0;

		for (int i = 0; i < 7; i++)
		{
			if (i == 3 || i == 5)
			{
				q_desired_(i) = DyrosMath::cubic(play_time_, control_start_time_, control_start_time_ + duration, q_init_(i), target_q_(i), 0, 0);
				q_desired_dot_(i) = DyrosMath::cubicDot(play_time_, control_start_time_, control_start_time_ + duration, q_init_(i), target_q_(i), 0, 0, hz_);
			}
			else
			{
				q_desired_(i) = 0;
				q_desired_dot_(i) = 0;
			}
		}
		
		Kp = 400 * EYE(7);
		Kv = 40 * EYE(7);
		torque_desired_ = m_ * (Kp * (q_desired_ - q_) + Kv * (q_desired_dot_ - qdot_)) + g_;;
	}
	else if (control_mode_ == "hw5-1-1")
	{
		Vector7d target_q_;
		target_q_ << 0.0, 0.0, 0.0, -25. * DEG2RAD, 0.0, M_PI / 2, 0.0;
		double duration = 3.0;

		Matrix6d lambda_;
		lambda_ = (j_ * m_.inverse() * j_.transpose()).inverse();
		Matrix<double, 7, 6> j_bar_;
		j_bar_ = m_.inverse() * j_.transpose() * lambda_;
		
		/*x_target_(0) = x_init_(0);
		x_target_(1) = x_init_(1) + 0.02;
		x_target_(2) = x_init_(2);*/

		if (play_time_ - control_start_time_ < 1)
		{
			x_target_(0) = x_init_(0);
			x_target_(1) = x_init_(1);
			x_target_(2) = x_init_(2);
		}
		else
		{
			x_target_(0) = x_init_(0);
			x_target_(1) = x_init_(1) + 0.02;
			x_target_(2) = x_init_(2);
		}

		Vector3d x_desired_dot_;
		x_desired_dot_ << 0, 0, 0;

		Kp_val = 400;
		Kv_val = 40;

		Vector3d F_star_;
		F_star_ = Kp_val * EYE(3) * (x_target_ - x_) + Kv_val * EYE(3) * (x_desired_dot_ - x_dot_.head(3));

		Vector3d M_star_;
		M_star_ = -Kp_val * EYE(3) * DyrosMath::getPhi(rotation_, rotation_init_) - Kv_val * x_dot_.tail<3>();

		Vector6d F_zero_star_;
		F_zero_star_ << F_star_, M_star_;

		Vector7d torque_zero_;
		torque_zero_ = m_ * (Kp_val * EYE(7) * (q_init_ - q_) - Kv_val * EYE(7) * qdot_);

		torque_desired_ = j_.transpose() * lambda_ * F_zero_star_ + (EYE(7) - j_.transpose() * j_bar_.transpose()) * torque_zero_ + g_;
		
		//fprintf(fp1, "%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t \n", play_time_ - control_start_time_, q_desired_(0), q_desired_(1), q_desired_(2), q_desired_(3), q_desired_(4), q_desired_(5), q_desired_(6), q_(0), q_(1), q_(2), q_(3), q_(4), q_(5), q_(6));
	}
	else if (control_mode_ == "hw5-1-2")
	{
		double duration = 3.0;

		Matrix6d lambda_;
		lambda_ = (j_ * m_.inverse() * j_.transpose()).inverse();
		Matrix<double, 7, 6> j_bar_;
		j_bar_ = m_.inverse() * j_.transpose() * lambda_;

		x_target_(0) = x_init_(0);
		x_target_(1) = x_init_(1) + 0.1;
		x_target_(2) = x_init_(2);

		Vector3d x_desired_dot_;
		x_desired_dot_ << 0, 0, 0;
		
		for (int i = 0; i < 3; i++) {
			x_desired_(i) = DyrosMath::cubic(play_time_, control_start_time_, control_start_time_ + duration, x_init_(i), x_target_(i), 0, 0);
			x_desired_dot_(i) = DyrosMath::cubicDot(play_time_, control_start_time_, control_start_time_ + duration, x_init_(i), x_target_(i), 0, 0, hz_);
		}

		Kp_val = 400;
		Kv_val = 40;

		Vector3d F_star_;
		F_star_ = Kp_val * EYE(3) * (x_desired_.head(3) - x_) + Kv_val * EYE(3) * (x_desired_dot_ - x_dot_.head(3));

		Vector3d M_star_;
		M_star_ = -Kp_val * EYE(3) * DyrosMath::getPhi(rotation_, rotation_init_) - Kv_val * x_dot_.tail<3>();

		Vector6d F_zero_star_;
		F_zero_star_ << F_star_, M_star_;

		Vector7d torque_zero_;
		torque_zero_ = m_ * (Kp_val * EYE(7) * (q_init_ - q_) - Kv_val * EYE(7) * qdot_);

		torque_desired_ = j_.transpose() * lambda_ * F_zero_star_ + (EYE(7) - j_.transpose() * j_bar_.transpose()) * torque_zero_ + g_;

		//fprintf(fp1, "%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t \n", play_time_ - control_start_time_, q_desired_(0), q_desired_(1), q_desired_(2), q_desired_(3), q_desired_(4), q_desired_(5), q_desired_(6), q_(0), q_(1), q_(2), q_(3), q_(4), q_(5), q_(6));
	}
	else if (control_mode_ == "hw5-1-2")
	{
		double duration = 3.0;

		Matrix6d lambda_;
		lambda_ = (j_ * m_.inverse() * j_.transpose()).inverse();
		Matrix<double, 7, 6> j_bar_;
		j_bar_ = m_.inverse() * j_.transpose() * lambda_;

		x_target_(0) = x_init_(0);
		x_target_(1) = x_init_(1) + 0.1;
		x_target_(2) = x_init_(2);

		Vector3d x_desired_dot_;
		x_desired_dot_ << 0, 0, 0;
		
		for (int i = 0; i < 3; i++) {
			x_desired_(i) = DyrosMath::cubic(play_time_, control_start_time_, control_start_time_ + duration, x_init_(i), x_target_(i), 0, 0);
			x_desired_dot_(i) = DyrosMath::cubicDot(play_time_, control_start_time_, control_start_time_ + duration, x_init_(i), x_target_(i), 0, 0, hz_);
		}

		Kp_val = 400;
		Kv_val = 40;

		Vector3d F_star_;
		F_star_ = Kp_val * EYE(3) * (x_desired_.head(3) - x_) + Kv_val * EYE(3) * (x_desired_dot_ - x_dot_.head(3));

		Vector3d M_star_;
		M_star_ = -Kp_val * EYE(3) * DyrosMath::getPhi(rotation_, rotation_init_) - Kv_val * x_dot_.tail<3>();

		Vector6d F_zero_star_;
		F_zero_star_ << F_star_, M_star_;

		Vector7d torque_zero_;
		torque_zero_ = m_ * (Kp_val * EYE(7) * (q_init_ - q_) - Kv_val * EYE(7) * qdot_);

		torque_desired_ = j_.transpose() * lambda_ * F_zero_star_ + (EYE(7) - j_.transpose() * j_bar_.transpose()) * torque_zero_ + g_;

		//fprintf(fp1, "%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t \n", play_time_ - control_start_time_, q_desired_(0), q_desired_(1), q_desired_(2), q_desired_(3), q_desired_(4), q_desired_(5), q_desired_(6), q_(0), q_(1), q_(2), q_(3), q_(4), q_(5), q_(6));
	}
	else if (control_mode_ == "hw5-2")
	{
		double duration = 3.0;

		Matrix6d lambda_;
		lambda_ = (j_ * m_.inverse() * j_.transpose()).inverse();
		Matrix<double, 7, 6> j_bar_;
		j_bar_ = m_.inverse() * j_.transpose() * lambda_;

		x_target_(0) = 0.3;
		x_target_(1) = -0.012;
		x_target_(2) = 0.52;

		Vector3d x_desired_dot_;
		x_desired_dot_ << 0, 0, 0;

		for (int i = 0; i < 3; i++) {
			x_desired_(i) = x_target_(i);
			x_desired_dot_(i) = 0;
		}

		Kp_val = 400.0;
		Kv_val = 40.0;

		Vector3d saturation_x_desired_dot;
		for (int i = 0; i < 3; i++) {
			saturation_x_desired_dot(i) = Kp_val / Kv_val * (x_desired_(i) - x_(i));
			if (abs(saturation_x_desired_dot(i)) > 0.3) {
				saturation_x_desired_dot(i) = 0.3 / abs(x_desired_(i) - x_(i)) * (x_desired_(i) - x_(i));
			}
		}

		Vector3d F_star_;
		F_star_ = Kv_val * EYE(3) * (saturation_x_desired_dot - x_dot_.head(3));

		Vector3d M_star_;
		M_star_ = -Kp_val * EYE(3) * DyrosMath::getPhi(rotation_, rotation_init_) - Kv_val * x_dot_.tail<3>();

		Vector6d F_zero_star_;
		F_zero_star_ << F_star_, M_star_;

		Vector7d torque_zero_;
		torque_zero_ = m_ * (Kp_val * EYE(7) * (q_init_ - q_) - Kv_val * EYE(7) * qdot_);

		torque_desired_ = j_.transpose() * lambda_ * F_zero_star_ + (EYE(7) - j_.transpose() * j_bar_.transpose()) * torque_zero_ + g_;

		//fprintf(fp1, "%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t \n", play_time_ - control_start_time_, q_desired_(0), q_desired_(1), q_desired_(2), q_desired_(3), q_desired_(4), q_desired_(5), q_desired_(6), q_(0), q_(1), q_(2), q_(3), q_(4), q_(5), q_(6));
	}
	else if (control_mode_ == "gravity")
	{
		torque_desired_ = g_;
	}
	else
	{
		Vector7d target_q_;
		target_q_ << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
		double duration = 1.0;

		for (int i = 0; i < 7; i++)
		{
			q_desired_(i) = DyrosMath::cubic(play_time_, control_start_time_, control_start_time_ + duration, q_init_(i), target_q_(i), 0, 0);
			q_desired_dot_(i) = DyrosMath::cubicDot(play_time_, control_start_time_, control_start_time_ + duration, q_init_(i), target_q_(i), 0, 0, hz_);
		}

		/*Kp = 300 * EYE(7);
		Kv = 1 * EYE(7);*/

		torque_desired_ = m_ * (Kp * (q_desired_ - q_) + Kv * (q_desired_dot_ - qdot_)) + g_;;
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

		cout << "play_time_ - control_start_time_ : ";
		cout << play_time_ - control_start_time_ << endl;
		cout << "q now    :\t";
		cout << std::fixed << std::setprecision(3) << q_.transpose() << endl;
		cout << "q desired:\t";
		cout << std::fixed << std::setprecision(3) << q_desired_.transpose() << endl;
		cout << "torque   :\t";
		cout << std::fixed << std::setprecision(3) << torque_.transpose() << endl;
		cout << "torque desired:\t";
		cout << std::fixed << std::setprecision(3) << torque_desired_.transpose() << endl;
		cout << "x_init_:\t";
		cout << std::fixed << std::setprecision(3) << x_init_.transpose() << endl;
		cout << "x_target_:\t";
		cout << std::fixed << std::setprecision(3) << x_target_.transpose() << endl;
		//cout << "m_:\n";
		//cout << std::fixed << std::setprecision(3) << m_ << endl;
		cout << "x        :\t";
		cout << x_.transpose() << endl;
		//cout << "x2        :\t";
		//cout << CalcBodyToBaseCoordinates(*model_, q_, body_id_[DOF - 4], com_position_[DOF - 4], false).transpose() << endl;
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

