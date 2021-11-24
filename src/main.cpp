
// Copyright (c) 2018 Seoul National University - DYROS
// Use of this source code is governed by the Apache-2.0 license, see LICENSE

// Native C
#include <sys/mman.h>	

#include <iostream>
#include <string>
#include <memory>
#include <thread>

#include <franka/exception.h>
#include <franka/robot.h>
#include <franka/model.h>

#include "controller.h"
#include "motion_generator.h"
#include "terminal_func.h"

using namespace std;

#define MODE(X,Y) case X: ac_ptr->setMode(Y); break;

bool running = true;

ArmController *ac_ptr;
std::mutex control_mode_mutex;

void inputCollector()
{
  while(running)
  {
    if (kbhit())
    {
      int key = getchar();
      control_mode_mutex.lock();
      switch (key)
      {
        // TODO: Implement with user input
        MODE('i', "joint_ctrl_init")
        MODE('h', "joint_ctrl_home")
        MODE('1', "simple_jacobian")
        MODE('2', "feedback_jacobian")
        MODE('3', "joint_ctrl_init")
        MODE('4', "CLIK")
        MODE('5', "Dynamics")
        MODE('6', "NullSpace")
        MODE('7', "VelocitySaturation")
      default:
        break;
      }
      control_mode_mutex.unlock();
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }
}

int main()
{
  if (mlockall(MCL_CURRENT | MCL_FUTURE))
    perror("mlockall failed:");
  
  franka::Robot robot("172.16.3.3");
  franka::Model model = robot.loadModel();

  robot.setCollisionBehavior(
      {{20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}}, {{20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}},
      {{10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0}}, {{10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0}},
      {{20.0, 20.0, 20.0, 20.0, 20.0, 20.0}}, {{20.0, 20.0, 20.0, 20.0, 20.0, 20.0}},
      {{10.0, 10.0, 10.0, 10.0, 10.0, 10.0}}, {{10.0, 10.0, 10.0, 10.0, 10.0, 10.0}});
  // robot.setJointImpedance({{3000, 3000, 3000, 2500, 2500, 2000, 2000}});
  robot.setJointImpedance({{2000, 2000, 2000, 1500, 1500, 1000, 1000}});
  // robot.setCartesianImpedance({{3000, 3000, 3000, 300, 300, 300}});

  // First move the robot to a suitable joint configuration
  std::array<double, 7> q_goal = {{0.0, 0.0, 0.0, -M_PI / 2., 0.0, M_PI / 2, 0}};
  auto q_goal_eigen = Eigen::Matrix<double, 7, 1>::Map(q_goal.data());

  MotionGenerator motion_generator(0.2, q_goal);
  std::cout << "WARNING: This program will move the robot! "
            << "Please make sure to have the user stop button at hand!" << std::endl
            << "Press Enter to continue..." << std::endl;
  std::cin.ignore();
  robot.control(motion_generator);
  std::cout << "Finished moving to initial joint configuration." << std::endl;

  usleep(1000000);

  const double hz = 1000.;
  const double period = 1./hz;
  ArmController ac(hz, model);
  double current_time = 0.0;
  ac_ptr = &ac;

  std::thread input_thread(inputCollector);


  bool is_first = true;
  // define callback for the position control loop
  std::function<franka::JointPositions(const franka::RobotState&, franka::Duration)>
      position_control_callback = [&](const franka::RobotState& robot_state,
                                        franka::Duration duration) -> franka::JointPositions {

    Eigen::Map<const Eigen::Matrix<double, 7, 1> > q(robot_state.q.data());
    Eigen::Map<const Eigen::Matrix<double, 7, 1> > dq(robot_state.dq.data());

    control_mode_mutex.lock();
    ac.readData(q,dq);
    if (is_first) {
      ac.readData(q_goal_eigen, Eigen::Vector7d::Zero());
      ac.initPosition(robot_state);
      is_first = false;
    }
    current_time += duration.toSec();
    
    ac.updateTime(duration.toSec());
    ac.compute();

    std::array<double, 7> position_d_array{};
    Eigen::Matrix<double, 7, 1>::Map(position_d_array.data()) = ac.getDesiredPosition();
    control_mode_mutex.unlock();
    franka::JointPositions output(position_d_array);
    return output;
  };

  
  // define callback for the torque control loop
  std::function<franka::Torques(const franka::RobotState&, franka::Duration)>
      torque_control_callback = [&](const franka::RobotState& robot_state,
                                        franka::Duration duration) -> franka::Torques {
    // get state variables
    std::array<double, 7> coriolis_array = model.coriolis(robot_state);

    // convert to Eigen
    Eigen::Map<const Eigen::Matrix<double, 7, 1> > coriolis(coriolis_array.data());
    Eigen::Map<const Eigen::Matrix<double, 7, 1> > q(robot_state.q.data());
    Eigen::Map<const Eigen::Matrix<double, 7, 1> > dq(robot_state.dq.data());

    control_mode_mutex.lock();
    ac.readData(q,dq);
    if (is_first) {
      ac.readData(q_goal_eigen, Eigen::Vector7d::Zero());
      ac.initPosition(robot_state);
      is_first = false;
    }
    current_time += period;
    
    ac.updateTime(duration.toSec());
    ac.compute();

    std::array<double, 7> tau_d_array{};
    Eigen::Matrix<double, 7, 1>::Map(&tau_d_array[0], 7) = ac.getDesiredTorque();
    control_mode_mutex.unlock();
    return tau_d_array;
  };

  try {
    // start real-time control loop
    std::cout << "WARNING: Collision thresholds are set to high values. "
              << "Make sure you have the user stop at hand!" << std::endl
              << "After starting try to push the robot and see how it reacts." << std::endl
              << "Press Enter to continue..." << std::endl;
    std::cin.ignore();

    // TODO: position mode : 
    // robot.control(position_control_callback);
    robot.control(position_control_callback, franka::ControllerMode::kJointImpedance, true, 10.0);

    // TODO: torque mode : 
    // robot.control(torque_control_callback);
  } catch (const franka::Exception& e) {
    std::cout << e.what() << std::endl;
    running = false;
    input_thread.join();
    return -1;
  }
  return 0;
}
