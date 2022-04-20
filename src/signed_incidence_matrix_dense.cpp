#include "signed_incidence_matrix_dense.h"

void signed_incidence_matrix_dense(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::MatrixXd & A)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  // Construct the sparse incidence matrix for a given edge network.
  //
  // Inputs: 
  //   n  number of vertices (#V)
  //   E  #E by 2 list of edge indices into rows of V
  // Outputs:
  //   A  #E by n signed incidence matrix
  A = Eigen::MatrixXd::Zero(E.rows(),n);
  int E_i, E_j;
  for(int i = 0; i < E.rows(); i++) {
    E_i = E.row(i).x();
    E_j = E.row(i).y();
    A(i, E_i) = 1;
    A(i, E_j) = -1;
  }

  // for(int i = 0; i < E.rows(); i++) {
  //   for(int j = 0; j < n; j++) {
  //     if(j == E.row(i).x()) A(i, j) = 1;
  //     else if(j == E.row(i).y()) A(i, j) = -1;
  //   }
  // }

  //////////////////////////////////////////////////////////////////////////////
}
