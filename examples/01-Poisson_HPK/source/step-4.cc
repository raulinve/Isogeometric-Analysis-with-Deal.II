/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 */


// @sect3{Include files}
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/convergence_table.h>  // IGA new

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>        // IGA
#include <deal.II/lac/sparse_direct.h>        // IGA
#include <deal.II/lac/precondition.h>         // IGA

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>             // IGA
#include <deal.II/fe/fe_bernstein.h>          // IGA
#include <deal.II/fe/mapping_fe_field.h>      // IGA

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>    // IGA

#include <fstream>
#include <iostream>
#include <list>                               // IGA
#include <string>                             // NEW
#include <filesystem>                         // NEW

#include "grid_generator.h"                   // IGA HEADER
#include "iga_handler.h"                      // IGA HEADER

using namespace dealii;

  bool original_problem_4 = false;	// Default: "false" [SHOULD BE REMOVED NEXT]

//====================================================
template <int dim>
class Solution : public Function<dim>
{
public:
  Solution () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  virtual Tensor<1, dim> gradient (const Point<dim>   &p,
                                   const unsigned int  component = 0) const;
};


/// This method implement the closed form solution of the analyzed problem.
template <int dim>
double Solution<dim>::value (const Point<dim>   &p,
                             const unsigned int) const
{
  return std::sin(2*p(0)*numbers::PI)*std::sin(3*p(1)*numbers::PI);
}


/// This method implement the gradient of the closed form solution of the analyzed problem. <br>
/// @note This method is used only to compute the H1 norm.
template <int dim>
Tensor<1, dim> Solution<dim>::gradient (const Point<dim>   &p,
                                const unsigned int) const
{
  // Rank 1 Tensor = Vector with (dim)-components;
  Tensor<1, dim> return_gradient;
  return_gradient[0] = (2*numbers::PI*std::cos(2*p(0)*numbers::PI)*std::sin(3*p(1)*numbers::PI));
  return_gradient[1] = (3*numbers::PI*std::sin(2*p(0)*numbers::PI)*std::cos(3*p(1)*numbers::PI));	
  return return_gradient;
}

//====================================================
/**
* @class   BoundaryValues
* @brief   Class used to define the solution at the boundaries.
*
* This class implement the function that the solution must have at the boundaries.
*/
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};


/// This method implement the function that the solution must have at the boundaries.
template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
  if (original_problem_4) {
    return p.square();
  }
  else {
    return std::sin(2*p(0)*numbers::PI)*std::sin(3*p(1)*numbers::PI);
  }
}


//====================================================
/**
* @class   RightHandSide
* @brief   Class used to define the force vector.
* 
* This class implement the function that acts as forcing term.
*/
template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};


/// This method implement the function that acts as forcing term.
template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int /*component*/) const
{
  if (original_problem_4) {
    double return_value = 0.0;
    for (unsigned int i = 0; i < dim; ++i)
      return_value += 4.0 * std::pow(p(i), 4.0);
    return return_value;
  }
  else {
    return 13*std::pow(numbers::PI,2)*std::sin(2*p(0)*numbers::PI)*std::sin(3*p(1)*numbers::PI);
  }
}





//====================================================
//====================================================
// @sect3{The <code>Laplace</code> class template}

template <int dim>
class Laplace
{
public:
  Laplace(IgaHandler<dim,dim>  & iga_handler,
          ConvergenceTable     & convergence_table);
  ~Laplace ();

  void run(unsigned int cycle);

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results(unsigned int cycle) const;
  void process_solution (const unsigned int cycle);     // new
  void print_table (const unsigned int cycle, unsigned int n_cycle);          // new

  IgaHandler<dim,dim>  &iga_handler;                    // IGA HEADER

  unsigned int         degree;                          // IGA
  unsigned int         cg_iter;                         // IGA  [ ?? used only in solve ]

  Triangulation<dim>   &triangulation;
  FE_Bernstein<dim>    &fe;                             // IGA MODIFICATION
  DoFHandler<dim>      &dof_handler;

    MappingFEField<dim>  *mappingfe;                      // IGA
    IndexSet             active_set;

    AffineConstraints<double>   bspline_constraints;      // IGA NAME
  //SparsityPattern      sparsity_pattern;
    SparsityPattern      sparsity_bspline;                // IGA [ NEW from setup_system]
  //SparseMatrix<double> system_matrix;
    SparseMatrix<double> bspline_system_matrix;           // IGA NAME

  //Vector<double> solution;
  //Vector<double> system_rhs;
    Vector<double>       bspline_solution;
    Vector<double>       bspline_system_rhs;

  ConvergenceTable     & convergence_table;  // IGA

};


// @sect4{Laplace::Laplace}

template <int dim>
Laplace<dim>::Laplace(IgaHandler<dim,dim>  & iga_handler,
                      ConvergenceTable     & convergence_table)
    : 
    iga_handler       (iga_handler),
    degree            (iga_handler.degree),
    triangulation     (iga_handler.tria),
    fe                (iga_handler.fe),
    dof_handler       (iga_handler.dh),
    mappingfe         (iga_handler.map_fe),
    convergence_table (convergence_table)
{
    std::cout << " - Initialization of the problem." << std::endl;
}



template <int dim>
Laplace<dim>::~Laplace()
{
  bspline_system_matrix.clear();
  if (mappingfe)
    delete mappingfe;
}


// @sect4{Laplace::make_grid}

template <int dim>
void Laplace<dim>::make_grid()
{
  //GridGenerator::hyper_cube(triangulation, -1, 1);
  //triangulation.refine_global(4);

  std::cout << "   Number of active cells: " 
            << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: " << triangulation.n_cells()
            << std::endl;
}

// @sect4{Laplace::setup_system}

// This function looks exactly like in the previous example, although it
// performs actions that in their details are quite different if
// <code>dim</code> happens to be 3. The only significant difference from a
// user's perspective is the number of cells resulting, which is much higher
// in three than in two space dimensions!
template <int dim>
void Laplace<dim>::setup_system()
{
    std::cout << " - Setup the system" << std::endl;

    dof_handler.distribute_dofs(fe);

    //active_set.set_size (iga_handler.n_bspline);
    //active_set_vector.reinit(iga_handler.n_bspline);

    std::cout << "   Number of degrees of freedom: " 
              << dof_handler.n_dofs()
              << std::endl
              << "   (Number of degrees of freedom IGA: "
              << iga_handler.n_bspline << ")"
              << std::endl;

  //DynamicSparsityPattern  dsp(dof_handler.n_dofs());
  //DoFTools::make_sparsity_pattern(dof_handler, dsp);
  //sparsity_pattern.copy_from(dsp);
  //system_matrix.reinit(sparsity_pattern);
  //solution.reinit(dof_handler.n_dofs());
  //system_rhs.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern  bspline_sp(iga_handler.n_bspline);
    iga_handler.make_sparsity_pattern (bspline_sp);
    sparsity_bspline.copy_from(bspline_sp);
    bspline_system_matrix.reinit(sparsity_bspline);
    bspline_solution.reinit (iga_handler.n_bspline);
    bspline_system_rhs.reinit (iga_handler.n_bspline);


//----------------------
// Boundary values
    QGauss<dim-1>  boundary_quad(fe.degree+2);

    std::map<types::global_dof_index,double> boundary_values;

    //typename FunctionMap<dim>::type  dirichlet_boundary;                            // <- [ DEPRECATED ]
    typename std::map<types::boundary_id, const Function<dim>*> dirichlet_boundary;   // <- [ NEW ]

    BoundaryValues<dim> boundary_funct;
    dirichlet_boundary[0] = &boundary_funct;

    iga_handler.project_boundary_values(dirichlet_boundary,
                                        boundary_quad,
                                        bspline_constraints);
//----------------------


}


// @sect4{Laplace::assemble_system}


template <int dim>
void Laplace<dim>::assemble_system()
{
  std::cout << " > ASSEMBLING THE SYSTEM (wait) ... " << std::endl;

    bspline_system_matrix = 0;
    bspline_system_rhs    = 0;

  //QGauss<dim> quadrature_formula(fe.degree + 1);
  //RightHandSide<dim> right_hand_side;
    const QGauss<dim>         quadrature_formula(fe.degree+1);
    const RightHandSide<dim>  right_hand_side;


  FEValues<dim> fe_values(*mappingfe, fe,
                          quadrature_formula,
                          update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Next, we again have to loop over all cells and assemble local
  // contributions.  Note, that a cell is a quadrilateral in two space
  // dimensions, but a hexahedron in 3D. In fact, the
  // <code>active_cell_iterator</code> data type is something different,
  // depending on the dimension we are in, but to the outside world they look
  // alike and you will probably never see a difference. In any case, the real
  // type is hidden by using `auto`:
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  //for (const auto &cell : dof_handler.active_cell_iterators())
    for (; cell!=endc; ++cell)
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      // Now we have to assemble the local matrix and right hand side.
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            //for (const unsigned int j : fe_values.dof_indices())
              for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_point) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_point) * // grad phi_j(x_q)
                 fe_values.JxW(q_point));           // dx

            const auto x_q = fe_values.quadrature_point(q_point);
            cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                            right_hand_side.value(x_q) *        // f(x_q)
                            fe_values.JxW(q_point));            // dx
          }
      // As a final remark to these loops: when we assemble the local
      // contributions into <code>cell_matrix(i,j)</code>, we have to multiply
      // the gradients of shape functions $i$ and $j$ at point number
      // q_index and
      // multiply it with the scalar weights JxW. This is what actually
      // happens: <code>fe_values.shape_grad(i,q_index)</code> returns a
      // <code>dim</code> dimensional vector, represented by a
      // <code>Tensor@<1,dim@></code> object, and the operator* that
      // multiplies it with the result of
      // <code>fe_values.shape_grad(j,q_index)</code> makes sure that the
      // <code>dim</code> components of the two vectors are properly
      // contracted, and the result is a scalar floating point number that
      // then is multiplied with the weights. Internally, this operator* makes
      // sure that this happens correctly for all <code>dim</code> components
      // of the vectors, whether <code>dim</code> be 2, 3, or any other space
      // dimension; from a user's perspective, this is not something worth
      // bothering with, however, making things a lot simpler if one wants to
      // write code dimension independently.

      // With the local systems assembled, the transfer into the global matrix
      // and right hand side is done exactly as before, but here we have again
      // merged some loops for efficiency:

      /*cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values.dof_indices())
        {
          for (const unsigned int j : fe_values.dof_indices())
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }*/

        iga_handler.distribute_local_to_global(cell_matrix,
                                               cell_rhs,
                                               cell,
                                               bspline_system_matrix,
                                               bspline_system_rhs,
                                               bspline_constraints);

    }

  // As the final step in this function, we wanted to have non-homogeneous
  // boundary values in this example, unlike the one before. This is a simple
  // task, we only have to replace the Functions::ZeroFunction used there by an
  // object of the class which describes the boundary values we would like to
  // use (i.e. the <code>BoundaryValues</code> class declared above):
  //
  // The function VectorTools::interpolate_boundary_values() will only work
  // on faces that have been marked with boundary indicator 0 (because that's
  // what we say the function should work on with the second argument below).
  // If there are faces with boundary id other than 0, then the function
  // interpolate_boundary_values will do nothing on these faces. For
  // the Laplace equation doing nothing is equivalent to assuming that
  // on those parts of the boundary a zero Neumann boundary condition holds.
  
/*std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           BoundaryValues<dim>(),
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);*/
//----------------------
// Boundary values
/*std::cout << " E" << std::endl;
    QGauss<dim-1>  boundary_quad(fe.degree+2);

    std::map<types::global_dof_index,double> boundary_values;

    //typename FunctionMap<dim>::type  dirichlet_boundary;                            // <- [ DEPRECATED ]
    typename std::map<types::boundary_id, const Function<dim>*> dirichlet_boundary;   // <- [ NEW ]

    BoundaryValues<dim> boundary_funct;
    dirichlet_boundary[0] = &boundary_funct;
std::cout << " F" << std::endl;
    iga_handler.project_boundary_values(dirichlet_boundary,
                                        boundary_quad,
                                        bspline_constraints);
std::cout << " G" << std::endl;*/
//----------------------

}


// @sect4{Laplace::solve}

// Solving the linear system of equations is something that looks almost
// identical in most programs. In particular, it is dimension independent, so
// this function is copied verbatim from the previous example.
template <int dim>
void Laplace<dim>::solve()
{
  std::cout << " > SOLVING THE SYSTEM (wait) ... " << std::endl;

  SolverControl             solver_control(100000, 1e-14);             // default: (1000, 1e-12)
  SolverCG<Vector<double>>  solver(solver_control);

  solver.solve(bspline_system_matrix, 
               bspline_solution, 
               bspline_system_rhs, 
               PreconditionIdentity());

  // We have made one addition, though: since we suppress output from the
  // linear solvers, we have to print the number of iterations by hand.
  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;

  cg_iter = solver_control.last_step();    // [ ?? ]
}


// @sect4{Laplace::output_results}

// This function also does what the respective one did in step-3. No changes
// here for dimension independence either.
//
// Since the program will run both 2d and 3d versions of the Laplace solver,
// we use the dimension in the filename to generate distinct filenames for
// each run (in a better program, one would check whether <code>dim</code> can
// have other values than 2 or 3, but we neglect this here for the sake of
// brevity).
template <int dim>
void Laplace<dim>::output_results(unsigned int cycle) const
{
  std::cout << " - Saving/overwriting the .vtk result files." << std::endl;

  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);

    Vector<double> bspline_sol_dh(dof_handler.n_dofs());
    Vector<double> bspline_solution_vector(bspline_solution);
    iga_handler.transform_vector_into_fe_space(
      bspline_sol_dh,
      bspline_solution_vector);

  //data_out.add_data_vector(solution, "solution");
    data_out.add_data_vector (bspline_sol_dh, "solution");

  data_out.build_patches();

  std::string filename = "solution-2d-";
  filename += ('0' + cycle);
  filename += ".vtk";
  std::string relpath = "RESULTS/" + filename;
  std::ofstream output(relpath);
  data_out.write_vtk(output);
}




/// Main problem class: This method process the solution in order to compute the L2 and the H1 norm.  <br>
/*! In addition it prepares the convergence table.
*/
template <int dim>
void Laplace<dim>::process_solution(const unsigned int cycle)
{
  std::cout << " - Computing solution errors." << std::endl;
  Vector<float> difference_per_cell (triangulation.n_active_cells());
  //Quadrature<dim>      error_quad    = QGauss<dim>(degree+3);
  Quadrature<dim>      error_quad    = QGauss<dim>(2*fe.degree+1);

    Vector<double> bspline_sol_dh(dof_handler.n_dofs());
    Vector<double> bspline_solution_vector(bspline_solution);

    iga_handler.transform_vector_into_fe_space(
      bspline_sol_dh,
      bspline_solution_vector);

  // Evaluate L2-norm error:
  VectorTools::integrate_difference (dof_handler,
                                     bspline_sol_dh,
                                     Solution<dim>(),
                                     difference_per_cell,
                                     error_quad,
                                     VectorTools::L2_norm);
  /*const double L2_error = VectorTools::compute_global_error(triangulation, 
  							    difference_per_cell,
  							    VectorTools::L2_norm);*/
    const double L2_error = difference_per_cell.l2_norm();

  // Evaluate H1-norm error:                                          // IGA only
  VectorTools::integrate_difference (dof_handler,
                                     bspline_sol_dh,
                                     Solution<dim>(),
                                     difference_per_cell,
                                     error_quad,
                                     VectorTools::H1_seminorm); //VectorTools::H1_norm);
  /*const double H1_error = VectorTools::compute_global_error(triangulation, 
  							    difference_per_cell,
  							    VectorTools::H1_norm);*/
    const double H1_error = difference_per_cell.l2_norm();

    const QTrapez<1>     q_trapez;
    const QIterated<dim> q_iterated (q_trapez, 5);
    VectorTools::integrate_difference (dof_handler, bspline_sol_dh,
                                       Solution<dim>(),
                                       difference_per_cell,
                                       q_iterated,
                                       VectorTools::Linfty_norm);

    const double Linfty_error = difference_per_cell.linfty_norm();

    const unsigned int n_active_cells=triangulation.n_active_cells();
    const unsigned int n_dofs=dof_handler.n_dofs();

    convergence_table.add_value("cycle", cycle);
    convergence_table.add_value("cells", n_active_cells);
    convergence_table.add_value("dofs", n_dofs);
    convergence_table.add_value("bsplines", iga_handler.n_bspline);
    convergence_table.add_value("CG", cg_iter);
    convergence_table.add_value("memory_sys", (unsigned int)bspline_system_matrix.memory_consumption());
    convergence_table.add_value("L2", L2_error);
    convergence_table.add_value("H1", H1_error);
    convergence_table.add_value("Linfty", Linfty_error);

}




// @sect4{Laplace::run}

// This is the function which has the top-level control over everything. Apart
// from one line of additional output, it is the same as for the previous
// example.
template <int dim>
void Laplace<dim>::run(unsigned int cycle)
{
	std::cout << " - Problem RUN (" << cycle << ")" << std::endl;

    std::cout << " - Checking/creating the output directory \"RESULTS\". " << std::endl;
    std::string OutputFolderName = "RESULTS";
    std::filesystem::create_directories(OutputFolderName);

  //std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;

  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results(cycle);

  process_solution (cycle);
}


//====================================================
/**
* Global function
*
* This function sets up the convergence table and print it on screen, on a .txt file and on a .tex file (if needed).
*/
void print_table (dealii::ConvergenceTable & convergence_table,
                  bool               h_refinement,
                  bool               p_refinement,
                  bool               k_refinement)
{
std::cout << "\n\n - Preparing the output converging table." << std::endl;

  convergence_table.set_precision("L2", 3);
  convergence_table.set_precision("H1", 3);
  convergence_table.set_precision("Linfty", 3);
  convergence_table.set_scientific("L2", true);
  convergence_table.set_scientific("H1", true);
  convergence_table.set_scientific("Linfty", true);
  convergence_table.set_tex_caption("cells", "\\# cells");
  convergence_table.set_tex_caption("dofs", "\\# dofs");
  convergence_table.set_tex_caption("bsplines", "\\# B-splines");
  convergence_table.set_tex_caption("L2", "@f$L^2@f$-error");
  convergence_table.set_tex_caption("H1", "@f$H^1@f$-error");
  convergence_table.set_tex_caption("Linfty", "@f$L^\\infty@f$-error");
  convergence_table.set_tex_format("cells", "r");
  convergence_table.set_tex_format("dofs", "r");
  convergence_table.set_tex_format("CG", "r");
  convergence_table.set_tex_format("memory_sys", "r");

  // (1) - Print on screen:
  //if (cycle==n_cycle)
    {
      //std::cout << "\nALL CYCLES ERRORS : " << std::endl;
      
      std::cout << "\nALL CYCLES CONVERGENCE SUMMARY : " << std::endl;
      std::cout << "---------------------------------------------------------------------- " << std::endl;
      convergence_table.write_text(std::cout);
      std::cout << "---------------------------------------------------------------------- " << std::endl;
    }

  /*
  convergence_table.add_column_to_supercolumn("cycle", "n cells");
  convergence_table.add_column_to_supercolumn("cells", "n cells");
  std::vector<std::string> new_order;
  new_order.push_back("n cells");
  new_order.push_back("dofs");
  new_order.push_back("bsplines");
  new_order.push_back("CG");
  new_order.push_back("memory_sys");
  new_order.push_back("L2");
  new_order.push_back("H1");
  new_order.push_back("Linfty");
  convergence_table.set_column_order (new_order);

  std::cout << "\nALL CYCLES CONVERGENCE RATE : " << std::endl;
  std::cout << "---------------------------------------------------------------------- " << std::endl;
  convergence_table.write_text(std::cout);
  std::cout << "---------------------------------------------------------------------- " << std::endl;
  */

  // (2) - Save .txt file:
  std::string error_filename = "now_error_deg3_c0";
  if (p_refinement)
    error_filename += "_p_ref.txt";
  if (h_refinement)
    error_filename += "_h_ref.txt";
  if (k_refinement)
    error_filename += "_k_ref.txt";

  std::string relpath = "RESULTS/" + error_filename;
  std::ofstream error_table_file(relpath);
  //std::ofstream error_table_file(error_filename.c_str());
  convergence_table.write_text(error_table_file);

}


//====================================================
//====================================================
//====================================================
// @sect3{The <code>main</code> function}

int main (int argc, char *argv[])
{
  std::cout << " " << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= THE CODE IS RUNNING ======================" << std::endl;
  std::cout << "============================================================" << std::endl;

  // DEFAULT VALUES:
  //--------------------------------------------------------
    unsigned int n_cycle       = 5;       // 5  (6)
    unsigned int degree        = 0;       // 0
    bool         h_refinement  = true;    // true
    bool         p_refinement  = false;   // false
    bool         k_refinement  = false;   // false
    std::string  continuity    = "C1";    // "C1"
  //--------------------------------------------------------

	// Initial rekap (on-screen print) -------------
	std::cout << "\n  SELECTED OPTIONS : " << std::endl; 
    std::cout << "   Dimensions            : " << "2" << std::endl; 
	std::cout << "   Number of cycles      : " << n_cycle << std::endl; 
	std::cout << "   Polynomial degree     : " << degree << std::endl; 
	std::cout << "   h_ref | p_ref | k_ref : " << h_refinement << " | "
                                              << p_refinement << " | "
                                              << k_refinement << std::endl; 
	std::cout << "   Continuity            : " << continuity << std::endl;    
	std::cout << std::endl; 

	if (!p_refinement && !h_refinement && !k_refinement) {
		std::cout << "  ERROR: Select one type of refinement! " << std::endl; 
		return 1;
	}
    // ---------------------------------------------

try
{
  deallog.depth_console(0);
	  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  // Initializations before the loops:
  std::vector<unsigned int> subdivisions(1);
  ConvergenceTable          convergence_table;

  std::cout << "\n > STARTING THE CYCLES: ========================" << std::endl;
  for (unsigned int cycle=1; cycle<=n_cycle; ++cycle)
    {
	  std::cout << "\n\n\n CYCLE # " << cycle << " of " << n_cycle << "  ======================" << std::endl;

      // REFINEMENT:START ========================================
      if (h_refinement)  {
		  std::cout << " - Setup h-refinement" << std::endl;
          subdivisions[0] = std::pow(2, cycle);
          degree = 3;  }
      if (p_refinement)  {
		  std::cout << " - Setup p-refinement" << std::endl;
          subdivisions[0] = 100;
          degree = cycle;  }
      if (k_refinement)  {
		  std::cout << " - Setup k-refinement" << std::endl;
          subdivisions[0] = 100;
          degree = 4;  }

	  // KNOTS MATRIX:
      std::vector<std::vector<double> > knots(1, std::vector<double>(subdivisions[0]+1));
      for (unsigned int i=0; i<subdivisions[0]+1; ++i)
        knots[0][i] = -1+i*2./subdivisions[0];

      //MULTS MATRIX:
      std::vector<std::vector<unsigned int> > mults(1, std::vector<unsigned int>(subdivisions[0]+1, 1));
	  if(continuity == "C0") {
			std::cout << " - Setup C0 continuity" << std::endl;
			for (unsigned int i=0; i<subdivisions[0]; ++i)
			mults[0][i] = degree;  }
	  else if(continuity == "C1") {
			std::cout << " - Setup C1 continuity" << std::endl;
			for (unsigned int i=0; i<subdivisions[0]; ++i)
				mults[0][i] = degree-1;  }
	  else if(continuity == "C2") {
			std::cout << " - Setup C2 continuity" << std::endl;
			for (unsigned int i=0; i<subdivisions[0]; ++i)
				mults[0][i] = degree-2;  }
      if (k_refinement)
        for (unsigned int i=0; i<subdivisions[0]; ++i)
          mults[0][i] = cycle+1;

      // open knot vectors
      mults[0][0] = degree+1;
      mults[0][subdivisions[0]] = degree+1;
      mults.push_back(mults[0]);
      knots.push_back(knots[0]);
      // REFINEMENT:END ========================================


	  // SETUP COMPLETED - STARTING THE MAIN PROCEDURES:
	  std::cout << " - Assemble IgaHandler ...  " << std::endl;
	  //Timings::Chrono Timer;		// NEW (START TIMER)
      //Timer.start();				// NEW
        IgaHandler<2,2> iga_hd2(knots, mults, degree);
      //Timer.stop();					// NEW (STOP TIMER)
      //std::cout << "   (Time to assemble the IgaHandler: "
      //          << Timer.wallTime()/1000000 << " seconds)" << std::endl;

      //=======================================
      Laplace<2> laplace_problem_2d(iga_hd2, convergence_table);
      laplace_problem_2d.run (cycle);
      //=======================================

      // END: =========

	  std::cout << " - CYCLE # " << cycle << " successfully ended!" << std::endl;
    }

	print_table(convergence_table, h_refinement, p_refinement, k_refinement);

}


	catch (std::exception &exc)
	{
	  std::cerr << std::endl << std::endl
		        << "----------------------------------------------------"
		        << std::endl;
	  std::cerr << "Exception on processing: " << std::endl
		        << exc.what() << std::endl
		        << "Aborting!" << std::endl
		        << "----------------------------------------------------"
		        << std::endl;

	  return 1;
	}
	catch (...)
	{
	  std::cerr << std::endl << std::endl
		        << "----------------------------------------------------"
		        << std::endl;
	  std::cerr << "Unknown exception!" << std::endl
		        << "Aborting!" << std::endl
		        << "----------------------------------------------------"
		        << std::endl;
	  return 1;
	}

  std::cout << " \n" << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= CODE ENDED CORRECTLY =====================" << std::endl;
  std::cout << "============================================================\n" << std::endl;

  return 0;
}
