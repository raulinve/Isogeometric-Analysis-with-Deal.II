/*! ---------------------------------------------------------------------
* Copyright (C) 1999 - 2020 by the deal.II authors
* Copyright (C) 2021 by Raul Invernizzi
*
* This file has been modified from the example program step-5 of the
* deal.II library.
*
* The deal.II library is free software; you can use it, redistribute
* it, and/or modify it under the terms of the GNU Lesser General
* Public License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
* The full text of the license can be found in the file LICENSE at
* the top level of the deal.II distribution.
*
* ---------------------------------------------------------------------
*
* Final Author: Raul Invernizzi 2021
*/


// @sect3{Include files}

// Include files:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_ilu.h>     // IGA

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>       // step-5: read a triangulation from disk
#include <deal.II/grid/manifold_lib.h>  // step-5: object describing the boundary of circular domain

#include <deal.II/fe/fe.h>              // IGA
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>      // IGA
#include <deal.II/fe/fe_bernstein.h>    // IGA
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <filesystem>                   // IGA



using namespace dealii;

  bool original_problem_5 = true;	// Default: "true" [SHOULD BE REMOVED NEXT]

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
  if (original_problem_5) {
    return 0*p(0);
  }
  else {
    return 0*p(0);
  }
}


//====================================================
// @sect3{Working with nonconstant coefficients}

// In this example, we want to use a variable coefficient in the
// elliptic operator. Since we have a function which just depends on
// the point in space we can do things a bit more simply and use a plain
// function instead of inheriting from Function.

// This is the implementation of the coefficient function for a single
// point. We let it return 20 if the distance to the origin is less than 0.5,
// and 1 otherwise.
template <int dim>
double coefficient(const Point<dim> &p)
{
  if (original_problem_5) {
    if (p.square() < 0.5 * 0.5)
      return 20;
    else
      return 1;
  }
  else {
    return 1;
  }
}


//====================================================
// @sect3{The <code>Step5</code> class template}

// The main class is mostly as in other deal.II examples. The most visible
// change is that the function <code>make_grid</code> has been
// removed, since creating the grid is now done in the <code>run</code>
// function and the rest of its functionality is now in
// <code>setup_system</code>. Apart from this, everything is as before.
template <int dim>
class Step5
{
public:
  Step5(const std::string  fe_name,        // IGA: add params
        const std::string  quadrature_name,
        const unsigned int degree,
        const unsigned int n_cycles_low,
        const unsigned int n_cycles_up);
  ~Step5();                                // IGA

  void run();

private:
  void setup_system();
  void assemble_system();
  void solve();
  void output_results(const unsigned int cycle) const;

  const std::string    fe_name;            // IGA
  const std::string    quadrature_name;    // IGA
  const unsigned int   degree;             // IGA
  const unsigned int   n_cycles_low;       // IGA
  const unsigned int   n_cycles_up;        // IGA
  unsigned int         cg_iter;            // IGA

  Triangulation<dim>   triangulation;
  //FE_Q<dim>            fe;
  FiniteElement<dim>   *fe;                // IGA
  DoFHandler<dim>      dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;

  Quadrature<dim>      matrix_quad;        // IGA
  Quadrature<dim>      error_quad;         // IGA
  Quadrature<dim-1>    boundary_quad;      // IGA
};




// @sect3{The <code>Step5</code> class implementation}

// @sect4{Step5::Step5}

// This function is as before.
template <int dim>
Step5<dim>::Step5(const std::string  fe_name,
                  const std::string  quadrature_name,
                  const unsigned int degree,
                  const unsigned int n_cycles_low,
                  const unsigned int n_cycles_up)
  : 
  fe_name         (fe_name),
  quadrature_name (quadrature_name),
  degree          (degree),
  n_cycles_low    (n_cycles_low),
  n_cycles_up     (n_cycles_up),
  fe              (NULL),
  dof_handler     (triangulation)
{
  // IGA: NEW IMPLEMENTATION

  if (quadrature_name == "legendre")       // [DEFAULT]
    {
      matrix_quad   = QGauss<dim>(degree+1);
      boundary_quad = QGauss<dim-1>(degree+2);
      error_quad    = QGauss<dim>(degree+3);
    }
  else if (quadrature_name == "lobatto")
    {
      matrix_quad   = QGaussLobatto<dim>(degree+2);
      boundary_quad = QGaussLobatto<dim-1>(degree+3);
      error_quad    = QGaussLobatto<dim>(degree+4);
    }
  else
    AssertThrow(false, ExcMessage("Quadrature not Implemented"));

  if      (fe_name == "bernstein")
    fe = new FE_Bernstein<dim>(degree);
  else if (fe_name == "lagrange")          // [DEFAULT]
    fe = new FE_Q<dim>(degree);
  else if (fe_name == "lobatto")
    fe = new FE_Q<dim>(QGaussLobatto<1>(degree+1));
  else
    AssertThrow(false, ExcMessage("FE not supported"));

}

/// Main problem class: Destructor
template <int dim>
Step5<dim>::~Step5 ()
{
  dof_handler.clear();
  if (fe)
    delete fe;
}


// @sect4{Step5::setup_system}

// This is the function <code>make_grid</code> from the previous
// example, minus the generation of the grid. Everything else is unchanged:
template <int dim>
void Step5<dim>::setup_system()
{
  dof_handler.distribute_dofs(*fe);        // IGA: "*"

  std::cout << "   Number of degrees of freedom: " 
            << dof_handler.n_dofs()
            << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}



// @sect4{Step5::assemble_system}

// The first parts of the function are completely unchanged from other examples.
template <int dim>
void Step5<dim>::assemble_system()
{
  std::cout << "   > ASSEMBLING THE SYSTEM (wait) ... " << std::endl;

  //QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(*fe, matrix_quad,                 // IGA: *, matrix_quad
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe->dofs_per_cell;     // IGA: "->"
  const unsigned int n_q_points    = matrix_quad.size();    // IGA

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Next is the typical loop over all cells to compute local contributions
  // and then to transfer them into the global matrix and vector. The only
  // change in this part, compared to step-4, is that we will use the
  // <code>coefficient</code> function defined above to compute the
  // coefficient value at each quadrature point.

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell!=endc; ++cell)
  //for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0.;
      cell_rhs    = 0.;
      fe_values.reinit(cell);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      //for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          const double current_coefficient =
            coefficient<dim>(fe_values.quadrature_point(q_point));
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          //for (const unsigned int i : fe_values.dof_indices())
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
              //for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) +=
                  (current_coefficient *                            // a(x_q)
                   fe_values.shape_grad(i, q_point) *               // grad phi_i(x_q)
                   fe_values.shape_grad(j, q_point) *               // grad phi_j(x_q)
                   fe_values.JxW(q_point));                         // dx

              cell_rhs(i) += (fe_values.shape_value(i, q_point) *   // phi_i(x_q)
                              1.0 *                                 // f(x_q)
                              fe_values.JxW(q_point));              // dx
            }
        }

      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      //for (const unsigned int i : fe_values.dof_indices())
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          //for (const unsigned int j : fe_values.dof_indices())
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

  // With the matrix so built, we use zero boundary values again:
  std::map<types::global_dof_index, double> boundary_values;

  using type=std::map<types::boundary_id, const Function<dim>*>;
  type dirichlet_boundary;
  BoundaryValues<dim> boundary_funct;
  dirichlet_boundary[0] = &boundary_funct;

  //VectorTools::interpolate_boundary_values(...)
  VectorTools::project_boundary_values(dof_handler,             // IGA (different method)
                                       dirichlet_boundary,      // IGA: 0,
                                       boundary_quad,           // IGA: Functions::ZeroFunction<dim>(),
                                       boundary_values);

  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}


// @sect4{Step5::solve}

// The solution process again looks mostly like in the other examples. 
// However, we will now use here a preconditioned conjugate gradient
// algorithm. 
// We will use SSOR (symmetric successive overrelaxation), with a relaxation 
// factor of 1.2. 
// For this purpose, the SparseMatrix class has a function which does one SSOR step,
// and we need to package the address of this function together with the
// matrix on which it should act (which is the matrix to be inverted) and the
// relaxation factor into one object. The <code>PreconditionSSOR</code> class
// does this for us. (<code>PreconditionSSOR</code> class takes a template
// argument denoting the matrix type it is supposed to work on. The default
// value is <code>SparseMatrix@<double@></code>, which is exactly what we need
// here, so we simply stick with the default and do not specify anything in
// the angle brackets.)
//
// With this, the rest of the function is trivial: instead of the
// <code>PreconditionIdentity</code> object we have created before, we now use
// the preconditioner we have declared, and the CG solver will do the rest for
// us:
template <int dim>
void Step5<dim>::solve()
{
  SolverControl            solver_control(100000, 1e-14);    // default: (1000, 1e-12)
  SolverCG<Vector<double>> solver(solver_control);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  std::cout << "   Memory consumption " << system_matrix.memory_consumption()
            << " bytes" << std::endl;

  std::cout << "   > SOLVING THE SYSTEM (wait) ... " << std::endl;
  solver.solve(system_matrix, solution, system_rhs, preconditioner);

  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;

  cg_iter = solver_control.last_step();
}


// @sect4{Step5::output_results and setting output flags}

// Writing output to a file is mostly the same as for the other tutorials.
//
// The function writes the output in VTU format, a variation of the VTK format
// that requires less disk space because it compresses the data. Of course,
// there are many other formats supported by the DataOut class if you
// desire to use a program for visualization that doesn't understand
// VTK or VTU.
template <int dim>
void Step5<dim>::output_results(const unsigned int cycle) const
{
  std::cout << "   Saving/overwriting the .vtu result files." << std::endl;
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches();

  std::string filename = "solution-";
  filename += ('0' + cycle);
  filename += ".vtu";
  std::string relpath = "RESULTS/" + filename;

  //std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
  std::ofstream output (relpath);
  data_out.write_vtu(output);
}



// @sect4{Step5::run}

// The second to last thing in this program is the definition of the
// <code>run()</code> function. In contrast to the previous programs, we will
// compute on a sequence of meshes that after each iteration is globally
// refined. The function therefore consists of a loop over 6 cycles.
//
// In previous examples, we have already used some of the functions from the
// <code>GridGenerator</code> class. Here we would like to read a grid from a
// file where the cells are stored and which may originate from someone else,
// or may be the product of a mesh generator tool.
//
// In order to read a grid from a file, we generate an object of data type
// GridIn and associate the triangulation to it (i.e. we tell it to fill our
// triangulation object when we ask it to read the file). Then we open the
// respective file and initialize the triangulation with the data in the file.
template <int dim>
void Step5<dim>::run()
{
  // Initial rekap (on-screen print)
  std::cout << " \n                    polyn     quad  deg refinements" << std::endl;
  std::cout << "SELECTED OPTIONS : " << fe_name << " "
                                     << quadrature_name << " "
                                     << degree << " "
                                     << n_cycles_low << " "
                                     << n_cycles_up << "\n" << std::endl;
                                       
  std::cout << "Checking/creating the output directory \"RESULTS\". " << std::endl;
  std::string OutputFolderName = "RESULTS";
  std::filesystem::create_directories(OutputFolderName);

  std::cout << "Reading the grid file \"circle-grid.inp\"." << std::endl;
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream input_file("circle-grid.inp");

  // ExcInternalError is a globally defined exception, which may be thrown
  // whenever something is terribly wrong. Example:  <code>dim</code> is not equal to two.
  Assert(dim == 2, ExcInternalError());

  // If you like to use another input format, you have to use one of the other
  // <code>grid_in.read_xxx</code> function. (See the documentation of the
  // <code>GridIn</code> class to find out what input formats are presently
  // supported.)
  grid_in.read_ucd(input_file);

  // The grid in the file describes a circle. Therefore we have to use a
  // manifold object which tells the triangulation where to put new points on
  // the boundary when the grid is refined. 
  const SphericalManifold<dim> boundary;
  triangulation.set_all_manifold_ids_on_boundary(0);
  triangulation.set_manifold(0, boundary);

  std::cout << "Solving problem in " << dim << " space dimensions. \n\n" << std::endl;   

  // Perform the cycles:
  for (unsigned int cycle = n_cycles_low; cycle < n_cycles_up+1; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle != 0)
        triangulation.refine_global(1);

      // Now that we have a mesh for sure, we write some output and do all the
      // things that we have already seen in the previous examples.
      std::cout << "   Number of active cells: "  
                << triangulation.n_active_cells() 
                << std::endl                      
                << "   Total number of cells: "   
                << triangulation.n_cells()        
                << std::endl;

      setup_system();
      assemble_system();
      solve();
      output_results(cycle);
    }
}


// @sect3{The <code>main</code> function}

// The main function looks mostly like the one in the previous example, so we
// won't comment on it further:
int main (int argc, char **argv)
{
  std::cout << " " << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= THE CODE IS RUNNING ======================" << std::endl;
  std::cout << "============================================================" << std::endl;

  // DEFAULT VALUES:
  //--------------------------------------------------------
  char fe_name[]             = "lagrange";    // “lagrange”
  char quad_name[]           = "legendre";    // “legendre”
  unsigned int degree        = 1;             // 1
  unsigned int n_cycles_down = 0;             // 0
  unsigned int n_cycles_up   = 5;             // 5
  //--------------------------------------------------------
  // ./step-5 lagrange legendre 1 0 1
  // ./step-5 bernstein legendre 1 0 1

  char *tmp[3];
  tmp[0] = argv[0];
  tmp[1] = fe_name;
  tmp[2] = quad_name;

  if (argc == 1)    // Default parameters
    {
      argv = tmp;
    }
  else              // User-passed parameters
    {
      AssertThrow(argc == 6,
                  ExcMessage("Wrong number of arguments: 0 or 5"));
      AssertThrow(sscanf(argv[3], "%ud", &degree) == 1,
                  ExcMessage("Unrecognized argument 3"));
      AssertThrow(sscanf(argv[4], "%ud", &n_cycles_down) == 1,
                  ExcMessage("Unrecognized argument 4"));
      AssertThrow(sscanf(argv[5], "%ud", &n_cycles_up) == 1,
                  ExcMessage("Unrecognized argument 5"));
    }

  //=======================================
  Step5<2> laplace_problem_2d(argv[1], argv[2], degree, n_cycles_down, n_cycles_up);
  laplace_problem_2d.run();
  //=======================================


  std::cout << " \n" << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= CODE ENDED CORRECTLY =====================" << std::endl;
  std::cout << "============================================================\n" << std::endl;
  return 0;
}
