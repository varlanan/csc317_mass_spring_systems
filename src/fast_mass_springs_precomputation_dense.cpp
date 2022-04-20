#include "fast_mass_springs_precomputation_dense.h"
#include "signed_incidence_matrix_dense.h"
#include <Eigen/Dense>

bool fast_mass_springs_precomputation_dense(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::MatrixXd & M,
  Eigen::MatrixXd & A,
  Eigen::MatrixXd & C,
  Eigen::LLT<Eigen::MatrixXd> & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  // Precompute matrices and factorizations necessary for the "Fast Simulation of
  // Mass-Spring Systems" method.
  //
  // Inputs: 
  //   V  #V by 3 list of vertex positions
  //   E  #E by 2 list of edge indices into rows of V
  //   k  spring stiffness
  //   m  #V list of masses 
  //   b  #b list of "pinned"/fixed vertices as indices into rows of V
  //   delta_t  time step in seconds
  // Outputs:
  //   r  #E list of edge lengths
  //   M  #V by #V mass matrix
  //   A  #E by #V signed incidence matrix
  //   C  #b by #V selection matrix
  //   prefactorization  LLT prefactorization of energy's quadratic matrix
  Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(V.rows(),V.rows());
  /////////////////////////////////////////////////////////////////////////////
  
  signed_incidence_matrix_dense(V.rows(), E, A);

  r.resize(E.rows());
  Eigen::Vector3d E_i, E_j;
  for(int i = 0; i < E.rows(); i++) {
    E_i = V.row(E.row(i).x());
    E_j = V.row(E.row(i).y());
    r[i] = (E_i - E_j).norm();
  }

  M = Eigen::MatrixXd::Zero(V.rows(), V.rows());
  for(int i = 0; i < V.rows(); i++) {
    M(i, i) = m[i];
  }

  C = Eigen::MatrixXd::Zero(b.rows(), V.rows());
  for(int i = 0; i < b.rows(); i++) {
    C(i, b(i)) = 1;
  }

  double w = 1e10;
  Q = k * A.transpose() * A + pow(delta_t, -2) * M + w * C.transpose() * C;

  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}
