#include "fast_mass_springs_precomputation_sparse.h"
#include "signed_incidence_matrix_sparse.h"
#include <vector>

bool fast_mass_springs_precomputation_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::SparseMatrix<double>  & M,
  Eigen::SparseMatrix<double>  & A,
  Eigen::SparseMatrix<double>  & C,
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  //std::vector<Eigen::Triplet<double> > ijv;
  //for(int i = 0;i<n;i++) ijv.emplace_back(i,i,1);
  Eigen::SparseMatrix<double> Q(V.rows(), V.rows());
  //Q.setFromTriplets(ijv.begin(),ijv.end());
  /////////////////////////////////////////////////////////////////////////////

  signed_incidence_matrix_sparse(V.rows(), E, A);

  r.resize(E.rows());
  Eigen::Vector3d E_i, E_j;
  for(int i = 0; i < E.rows(); i++) {
    E_i = V.row(E.row(i).x());
    E_j = V.row(E.row(i).y());
    r[i] = (E_i - E_j).norm();
  }
  std::vector<Eigen::Triplet<double> > ijv_M;
  M.resize(V.rows(), V.rows());
  for(int i = 0; i < V.rows(); i++) {
    ijv_M.emplace_back(i, i, m[i]);
  }
  M.setFromTriplets(ijv_M.begin(), ijv_M.end());

  std::vector<Eigen::Triplet<double> > ijv_C;
  C.resize(b.rows(), V.rows());
  for(int i = 0; i < b.rows(); i++) {
    ijv_C.emplace_back(i, b(i), 1);
  }
  C.setFromTriplets(ijv_C.begin(), ijv_C.end());

  double w = 1e10;
  Q = k * A.transpose() * A + pow(delta_t, -2) * M + w * C.transpose() * C;

  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}
