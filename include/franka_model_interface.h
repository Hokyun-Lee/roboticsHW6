
#include <franka/robot.h> 
#include <franka/model.h>
#include <Eigen/Dense>

class FrankaModelInterface
{
public:
  FrankaModelInterface(franka::Model& model) : model_(model)
  {}

  // For mass/inertia matrix
  void setRobotState(franka::RobotState state)
  {
    state_ = state;
  }

  Eigen::Affine3d getTransform(
    franka::Frame frame,
    const Eigen::Matrix<double, 7, 1> &q)
  {
    state_q_ = q;
    return Eigen::Affine3d(Eigen::Matrix4d::Map(model_.pose(frame,state_).data()));
  }

  Eigen::Matrix<double, 6, 7> getJacobianMatrix(
    franka::Frame frame,
    const Eigen::Matrix<double, 7, 1> &q)
  {
    state_q_ = q;
    return Eigen::Matrix<double, 6, 7>::Map(model_.zeroJacobian(frame, state_).data());
  }

  Eigen::Matrix<double, 7, 7> getMassMatrix(
    const Eigen::Matrix<double, 7, 1> &q)
  {
    state_q_ = q;
    return Eigen::Matrix<double, 7, 7>::Map(model_.mass(state_).data());
  }

private:
  franka::Model& model_;
  franka::RobotState state_;
  Eigen::Map<Eigen::Matrix<double, 7, 1>> state_q_{state_.q.data()};

};