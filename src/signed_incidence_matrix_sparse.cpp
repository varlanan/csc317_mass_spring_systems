#include "signed_incidence_matrix_sparse.h"
#include <vector>

void signed_incidence_matrix_sparse(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::SparseMatrix<double>  & A)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  std::vector<Eigen::Triplet<double> > ijv;

  A.resize(E.rows(),n);

  int E_i, E_j;
  for(int i = 0; i < E.rows(); i++) {
    E_i = E.row(i).x();
    E_j = E.row(i).y();
    ijv.emplace_back(i, E_i, 1);
    ijv.emplace_back(i, E_j, -1);
  }

  A.setFromTriplets(ijv.begin(),ijv.end());
  //////////////////////////////////////////////////////////////////////////////
}
