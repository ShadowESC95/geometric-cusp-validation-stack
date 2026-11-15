// Fully specified deal.II 9.x 1D FEM cross-check.
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <fstream>
#include <cmath>
using namespace dealii;

class GaussianSource : public Function<1> {
public: double value(const Point<1>& p, const unsigned int=0) const override { return std::exp(-p[0]*p[0]); }
};

int main(){
  const double L=12.0, D=0.25, mass=0.7;
  Triangulation<1> tria; GridGenerator::subdivided_hyper_cube(tria,256,-L/2,L/2);
  FE_Q<1> fe(2); DoFHandler<1> dofs(tria); dofs.distribute_dofs(fe);
  QGauss<1> q(4); DynamicSparsityPattern dsp(dofs.n_dofs()); DoFTools::make_sparsity_pattern(dofs,dsp);
  SparsityPattern sp; sp.copy_from(dsp); SparseMatrix<double> A(sp); Vector<double> b(dofs.n_dofs());
  MatrixTools::create_mass_matrix(dofs,q,A);
  A *= mass*mass;
  SparseMatrix<double> K(sp); MatrixTools::create_laplace_matrix(dofs,q,K); K *= D; A += K;
  VectorTools::create_right_hand_side(dofs,q,GaussianSource(),b);
  std::map<types::global_dof_index,double> bc; VectorTools::interpolate_boundary_values(dofs,Functions::ZeroFunction<1>(),bc);
  Vector<double> S(dofs.n_dofs()); MatrixTools::apply_boundary_values(bc,A,S,b);
  SolverControl control(2000,1e-12); SolverCG<Vector<double>> solver(control); solver.solve(A,S,b,PreconditionSSOR<SparseMatrix<double>>());
  DataOut<1> data; data.attach_dof_handler(dofs); data.add_data_vector(S,"S_OB_proxy"); data.build_patches(); std::ofstream out("dealII_boundary_solution.vtu"); data.write_vtu(out);
  return 0;
}
