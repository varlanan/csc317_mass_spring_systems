#include "fast_mass_springs_step_sparse.h"
#include <igl/matlab_format.h>

void fast_mass_springs_step_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::SparseMatrix<double>  & M,
  const Eigen::SparseMatrix<double>  & A,
  const Eigen::SparseMatrix<double>  & C,
  const Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  igl::matlab_format();
  
  Eigen::MatrixXd y = Eigen::MatrixXd::Identity(V.rows(), 3);
  Eigen::MatrixXd d = Eigen::MatrixXd::Identity(E.rows(), 3);
  Eigen::MatrixXd penalty = Eigen::MatrixXd::Identity(V.rows(), 3);
  Eigen::Vector3d E_i, E_j;
  
  double w = 1e10;
  y = pow(delta_t, -2) * M * ( 2 * Ucur - Uprev ) + fext;
  penalty = w * C.transpose() * C * V;

  Unext = Ucur;
  for(int iter = 0;iter < 50;iter++)
  {
    d = Eigen::MatrixXd::Identity(E.rows(), 3);
    for(int i = 0; i < E.rows(); i++) {
      E_i = Unext.row(E.row(i).x());
      E_j = Unext.row(E.row(i).y());
      d.row(i) = r[i] * (E_i - E_j).normalized();
    }
    const Eigen::MatrixXd b = k * A.transpose() * d  + y + penalty;
    Unext = prefactorization.solve(b);
  }
  //////////////////////////////////////////////////////////////////////////////
}
