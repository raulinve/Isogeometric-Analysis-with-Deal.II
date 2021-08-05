/**
* @file         step-41.cc
* @brief        step-41 tutorial solved using isogeometric analysis.
* @detail       Isogeometric Analysis for the solution of an obstacle problem taken from the step-41 tutorial of the deal.II library.
* @author       Raul Invernizzi, Marco Tezzele, Nicola Cavallini, Luca Heltai.
*/
// @include      grid_generator.h iga_handler.h

/*! \mainpage MAIN PAGE
*
* Isogeometric Analysis classes for the deal.II library.
* 
* <hr>
* This code is used to introduce the Isogeometric Analysis (IGA) framework 
* into the standard deal.II Finite Element (FE) library.
* 
* This specific example named "step-41" is based on the Laplace equation 
* in 2d and deals with the question "what happens if a membrane is deflected 
* by some external force but is also constrained by an obstacle?" <br>
* 
* In other words, think of a elastic membrane clamped at the boundary 
* to a rectangular frame, its shape will be modified due to gravity acting on it.  <br>
* What happens now if there is an obstacle under the membrane that prevents 
* it from reaching its equilibrium position if gravity was the only existing force?  <br>
* 
* Please refere to the most recent version available on GitHub at the following link:  
* <a href="https://github.com/raulinve/Isogeometric-Analysis-with-Deal.II">
* Isogeometric-Analysis-with-Deal.II</a>. \n\n
*
* <hr>
* <b>Usage:</b>  
*
* The executable takes no arguments, to run it is sufficient to use the following command: <br> 
* <pre>
* <tt><b>    ./exe  </b></tt>
* </pre>
* 
* The default parameters are:: <br> 
* <pre>
* <tt>    Number of cycles       = 6      </tt>
* <tt>    Polynomial degree      = 0      </tt>
* <tt>    Type of refinement     = h-ref  </tt>
* <tt>    Polynomial continuity  = C1     </tt>
* </pre>
* 
* <hr>
* <b>Classes & Methods:</b>  
*
* int <b>main</b> (int argc, char **argv)<br> 
* 
* class <b>"Solution"</b>:<br>
* Function: @f$ val = {p(0)}^2 + {p(1)}^2 - 0.49  [if >0 ={val}^2]  [if <=0 =0] @f$ <br> 
* class Solution : public Function<dim>  
* double Solution<dim>::value (...)  
* Tensor<1, dim> Solution<dim>::gradient  
* 
* class <b>"BoundaryValues"</b>:<br> 
* Function: @f$ ( {p(0)}^2) + {p(1)}^2 - 0.49)^2 @f$ <br> 
* class BoundaryValues : public Function<dim>  
* double BoundaryValues<dim>::value  
* BoundaryValues::value (solution at the boundary) :  
* 
* class <b>"RightHandSide"</b>:<br> 
* Given: @f$ p.distance(Point) @f$  <br>
* @f$ if >  0.5 : -8\cdot (2\cdot {p(0)}^2 + 2\cdot {p(1)}^2 - 0.49) @f$  <br>
* @f$ if <= 0.5 : -8\cdot 0.49\cdot (1-{p(0)}^2 - {p(1)}^2 + 0.49)   @f$  <br>
* class RightHandSide : public Function<dim>  
* double RightHandSide<dim>::value  
* 
* class <b>"Obstacle"</b>:<br> 
* obstacle: @f$ =0 @f$  <br>
* class RightHandSide : public Function<dim>  
* double Obstacle<dim>::value
* 
* class <b>"Laplace"</b>:<br> 
* class Laplace  
* Laplace<dim>::Laplace  
* Laplace<dim>::~Laplace ()  
* void Laplace<dim>::make_grid ()  
* void Laplace<dim>::setup_system ()  
* void Laplace<dim>::assemble_system ()  
* void Laplace<dim>::solve ()  
* void Laplace<dim>::refine_grid()  
* void Laplace<dim>::output_results  
* void Laplace<dim>::process_solution  
* void Laplace<dim>::print_table  
* void Laplace<dim>::run ()  
* 
* <hr>

*/


/*! ---------------------------------------------------------------------
* Copyright (C) 1999 - 2015 by the deal.II authors
* Copyright (C) 2015 by Marco Tezzele, Nicola Cavallini, Luca Heltai
* Copyright (C) 2021 by Raul Invernizzi
*
* This file has been modified from the example program step-41 of the
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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/convergence_table.h>   // IGA

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
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
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>    // IGA

#include <fstream>
#include <iostream>
#include <list>                               // IGA
#include <string>                             // NEW
#include <filesystem>                         // NEW

#include "grid_generator.h"                   // IGA HEADER
#include "iga_handler.h"                      // IGA HEADER

#include "chrono.hpp"                         // NEW (timer)



namespace Step41
{
  using namespace dealii;

  bool original_problem_41 = false;		// Default: "false" [SHOULD BE REMOVED NEXT]

//====================================================
/**
* @class   Solution
* @brief   Class used to define the exact solution.
*
* This class implement the closed form solution of the analyzed problem.
* It contains a value method which implement the expression of the exact solution 
* and a gradient method used to implement the closed form gradient expression of the exact solution.
*/
  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    Solution () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                    const unsigned int  component = 0) const;
  };


/// This method implement the closed form solution of the analyzed problem.
  template <int dim>
  double Solution<dim>::value (const Point<dim> &p,
                               const unsigned int /*component*/) const
  {
    double value = p(0)*p(0) + p(1)*p(1) - 0.49;

    if ( value > 0 )
      return std::pow( value, 2);
    else
      return 0;
  }


/// This method implement the gradient of the closed form solution of the analyzed problem. <br>
/// @note This method is used only to compute the H1 norm.
  template <int dim>
  Tensor<1,dim> Solution<dim>::gradient (const Point<dim> &p,
                                         const unsigned int /*component*/) const
  {
    Tensor<1,dim> return_value;

    double value = p(0)*p(0) + p(1)*p(1) - 0.49;

    if ( value > 0)
      {
        return_value[0] = 4 * p(0) * value;
        return_value[1] = 4 * p(1) * value;
      }
    else
      {
        return_value[0] = 0;
        return_value[1] = 0;
      }
    return return_value;
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
                                    const unsigned int component) const
  {
     if (original_problem_41) {
		(void)component;
		AssertIndexRange(component, 1);
		return -10;
     } else {
		Assert (component == 0, ExcNotImplemented());

		if (p.distance(Point<dim>()) > 0.5)
		  return -8 * (2*p(0)*p(0) + 2*p(1)*p(1) - 0.49);
		else
		  return -8 * 0.49 * (1 - p(0)*p(0) - p(1)*p(1) + 0.49);
     }
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
                                     const unsigned int component) const
  {
     if (original_problem_41) {
        (void)component;
        AssertIndexRange(component, 1);
        return 0;
     } else {
        Assert (component == 0, ExcNotImplemented());
        return std::pow( p(0)*p(0) + p(1)*p(1) - 0.49, 2);
     }
  }


//====================================================
/**
* @class   Obstacle
* @brief   Class used to define the obstacle forcing the membrane.
* 
* This class implement the function characterizing the obstacle over which the membrane lies.
* 
*/
  template <int dim>
  class Obstacle : public Function<dim>
  {
  public:
    Obstacle () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };

  template <int dim>
  double Obstacle<dim>::value (const Point<dim> &p,
                               const unsigned int component) const
  {
     if (original_problem_41) {
        (void)component;
        Assert(component == 0, ExcIndexRange(component, 0, 1));
        if (p(0) < -0.5)
          return -0.2;
        else if (p(0) >= -0.5 && p(0) < 0.0)
          return -0.4;
        else if (p(0) >= 0.0 && p(0) < 0.5)
          return -0.6;
        else
          return -0.8;
     } else {
        Assert (component == 0, ExcNotImplemented());
        int val = p(0) * 0;    // <- [ USELESS - Just to avoid warnings of unused params ]
        return val*0;
     }
  }


//====================================================
/**
* @class   ObstacleProblem
* @brief   Main problem class
* 
* This class implement and collect every aspect of the defined problem.
*
* @note Differences with respect to the standard step-41. [TO BE UPDATED]
* 
*/
  template <int dim>
  class ObstacleProblem
  {
  public:
    ObstacleProblem (IgaHandler<dim,dim> &iga_handler,
                     ConvergenceTable &convergence_table);
    ~ObstacleProblem ();
    void run (unsigned int cycle);

  private:
    void make_grid ();
    void setup_system();
    void assemble_system ();
    void assemble_mass_matrix_diagonal (SparseMatrix<double> &mass_matrix);
    void update_solution_and_constraints ();
    void solve ();
    void output_results (const unsigned int iteration);
    void compute_error (unsigned int cycle);              // IGA


    IgaHandler<dim,dim>  &iga_handler;                    // IGA HEADER

    unsigned int         degree;                          // IGA
    unsigned int         cg_iter;                         // IGA  [ ?? used only in solve ]

    Triangulation<dim>   &triangulation;
    FE_Bernstein<dim>    &fe;                             // IGA MODIFICATION
    DoFHandler<dim>      &dof_handler;

    MappingFEField<dim>  *mappingfe;                      // IGA
    IndexSet             active_set;

    AffineConstraints<double>   bspline_constraints;      // IGA NAME

    SparseMatrix<double> bspline_system_matrix;           // IGA NAME
    SparseMatrix<double> bspline_complete_system_matrix;  // IGA NAME

    SparseMatrix<double> bspline_mass_matrix;             // IGA [ TAKE OUT from setup_system ]
    Vector<double>       mass_lumping;                    // IGA [ TAKE OUT from setup_system ]
    Vector<double>       diagonal_of_mass_matrix;
    Vector<double>       bspline_mass_lumping;            // IGA [ no ref here ?? ]

    Vector<double>       bspline_solution;
    Vector<double>       bspline_system_rhs;
    Vector<double>       bspline_complete_system_rhs;

    Vector<double>       bspline_contact_force;
    Vector<double>       bspline_lambda;                  // IGA [ TAKE OUT from update_solution_and_constraints ]

    SparsityPattern      sparsity_bspline;                // IGA [ NEW from setup_system]

    TrilinosWrappers::PreconditionAMG   precondition;     // IGA[ TAKE OUT from solve ]
    Vector<double>                 active_set_vector;     // IGA[ TAKE OUT from output_results ]

    ConvergenceTable              &convergence_table;     // IGA [ NEW Addition ]
  };

//-------------------------
/// Main problem class: Constructor
  template <int dim>
  ObstacleProblem<dim>::ObstacleProblem (IgaHandler<dim,dim> &iga_handler,
                                         ConvergenceTable &convergence_table)
    :
    iga_handler(iga_handler),
    degree(iga_handler.degree),
    triangulation(iga_handler.tria),
    fe (iga_handler.fe),
    dof_handler (iga_handler.dh),
    mappingfe(iga_handler.map_fe),
    convergence_table(convergence_table)
  {	
    std::cout << " - Initialization of the obstacle problem." << std::endl;
  }

//-------------------------
/// Main problem class: Destructor
  template <int dim>
  ObstacleProblem<dim>::~ObstacleProblem ()
  {
    bspline_system_matrix.clear();
    bspline_complete_system_matrix.clear();
    bspline_mass_matrix.clear();
    if (mappingfe)
      delete mappingfe;
  }

//-------------------------
/// Main problem class: Method used to generate the triangulation and to refine it.
  template <int dim>
  void ObstacleProblem<dim>::make_grid ()
  {
    //std::cout << " - Making the grid" << std::endl;
    std::cout << "   (Degree: "
              << degree << ")"
              << std::endl
              << "   (Number of active cells: "
              << triangulation.n_active_cells() << ")"
              << std::endl
              << "   (Total number of cells:  "
              << triangulation.n_cells() << ")"
              << std::endl;
  }

//-------------------------
/// Main problem class: Method used to setup the matrix system.
  template <int dim>
  void ObstacleProblem<dim>::setup_system ()
  {
    std::cout << " - Setup the system" << std::endl;

    dof_handler.distribute_dofs (fe);

    active_set.set_size (iga_handler.n_bspline);
    active_set_vector.reinit(iga_handler.n_bspline);

    std::cout << "   (Number of degrees of freedom:     "
              << dof_handler.n_dofs() << ")"
              << std::endl
              << "   (Number of degrees of freedom IGA: "
              << iga_handler.n_bspline << ")"
              << std::endl;

    DynamicSparsityPattern bspline_sp(iga_handler.n_bspline);
    iga_handler.make_sparsity_pattern (bspline_sp);

    sparsity_bspline.copy_from(bspline_sp);
    bspline_system_matrix.reinit(sparsity_bspline);
    bspline_complete_system_matrix.reinit(sparsity_bspline);

    bspline_solution.reinit (iga_handler.n_bspline);
    bspline_system_rhs.reinit (iga_handler.n_bspline);
    bspline_complete_system_rhs.reinit (iga_handler.n_bspline);

    bspline_contact_force.reinit (iga_handler.n_bspline);
    bspline_lambda.reinit (iga_handler.n_bspline);

    bspline_mass_matrix.reinit (sparsity_bspline);
    assemble_mass_matrix_diagonal(bspline_mass_matrix);
    diagonal_of_mass_matrix.reinit (iga_handler.n_bspline);

    // Instead of extracting the diagonal of the mass matrix, we perform a
    // mass lumping:
    mass_lumping.reinit(iga_handler.n_bspline);
    for (unsigned int i=0; i<mass_lumping.size(); ++i)
      mass_lumping[i] = 1;
    bspline_mass_matrix.vmult(diagonal_of_mass_matrix, mass_lumping);

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
  }

//-------------------------
/// Main problem class: Method used to assemble the main system.
  template <int dim>
  void ObstacleProblem<dim>::assemble_system ()
  {
  std::cout << "  |> ASSEMBLING THE SYSTEM (wait) ... " << std::endl;

    bspline_system_matrix = 0;
    bspline_system_rhs    = 0;

    const QGauss<dim>         quadrature_formula(fe.degree+1);
    const RightHandSide<dim>  right_hand_side;

    FEValues<dim>             fe_values (*mappingfe, fe, quadrature_formula,
                                         update_values   | update_gradients |
                                         update_quadrature_points |
                                         update_JxW_values);

    const unsigned int        dofs_per_cell = fe.dofs_per_cell;
    const unsigned int        n_q_points    = quadrature_formula.size();

    FullMatrix<double>        cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>            cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        cell_matrix = 0;
        cell_rhs = 0;

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                     fe_values.shape_grad (j, q_point) *
                                     fe_values.JxW (q_point));

              cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                              right_hand_side.value (fe_values.quadrature_point (q_point)) *
                              fe_values.JxW (q_point));
            }

        iga_handler.distribute_local_to_global(cell_matrix,
                                               cell_rhs,
                                               cell,
                                               bspline_system_matrix,
                                               bspline_system_rhs,
                                               bspline_constraints);
      }
  }

//-------------------------
/// Main problem class: Method used to assemble the mass matrix.
  template <int dim>
  void
  ObstacleProblem<dim>::
  assemble_mass_matrix_diagonal (SparseMatrix<double> &mass_matrix)
  {
    std::cout << " - Assembling the mass matrix ... " << std::endl;

    QIterated<dim> quadrature_formula(QTrapez<1>(),fe.degree);

    FEValues<dim>             fe_values (*mappingfe, fe,
                                         quadrature_formula,
                                         update_values   |
                                         update_JxW_values);

    const unsigned int        dofs_per_cell = fe.dofs_per_cell;
    const unsigned int        n_q_points    = quadrature_formula.size();

    FullMatrix<double>        cell_matrix (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        cell_matrix = 0;

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_value (i, q_point) *
                                   fe_values.shape_value (j, q_point) *
                                   fe_values.JxW (q_point));

        iga_handler.distribute_local_to_global(cell_matrix,
                                               cell,
                                               mass_matrix,
                                               bspline_constraints);
      }
  }

//-------------------------
/// Main problem class: Method used to update the solution.
  template <int dim>
  void
  ObstacleProblem<dim>::update_solution_and_constraints ()
  {
    std::cout << "  |- Updating active set..." << std::endl;

    const double penalty_parameter = 100.0;

    bspline_lambda = 0;
    bspline_complete_system_matrix.residual (bspline_lambda,
                                             bspline_solution,
                                             bspline_complete_system_rhs);

    for (unsigned int i=0; i<bspline_contact_force.size(); ++i)
      bspline_contact_force[i] = bspline_lambda[i]/diagonal_of_mass_matrix[i];
    bspline_contact_force *= -1;

    bspline_constraints.clear();
    active_set.clear();
    active_set_vector.reinit(0);
    active_set_vector.reinit(iga_handler.n_bspline);

    const Obstacle<dim> obstacle;
    std::vector<bool>   dof_touched (iga_handler.n_bspline, false);
    std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);

    std::vector<Point<dim> > bspline_support_points(iga_handler.n_bspline);

    iga_handler.map_dofs_to_support_points(bspline_support_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        local_dof_indices = iga_handler.iga_objects[cell].dof_indices;

        for (unsigned int v=0; v<fe.dofs_per_cell; ++v)
          {
            const unsigned int dof_index = local_dof_indices[v];

            if (dof_touched[dof_index] == false)
              dof_touched[dof_index] = true;
            else
              continue;


            const double obstacle_value = obstacle.value (bspline_support_points[dof_index]);
            const double solution_value = bspline_solution (dof_index);

            if (bspline_lambda(dof_index) +
                penalty_parameter *
                diagonal_of_mass_matrix(dof_index) *
                (solution_value - obstacle_value)
                <
                0)
              {
                active_set.add_index (dof_index);
                active_set_vector[dof_index] = 1;
                bspline_constraints.add_line (dof_index);
                bspline_constraints.set_inhomogeneity (dof_index, obstacle_value);

                bspline_solution (dof_index) = obstacle_value;
                bspline_lambda (dof_index) = 0;
              }
          }
      }

    std::cout << "  |  (Size of active set: " << active_set.n_elements() << ")"
              << std::endl;

    std::cout << "  |  (Residual of the non-contact part of the system: "
              << bspline_lambda.l2_norm() << ")"
              << std::endl;


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
  }

//-------------------------
/// Main problem class: This method is called in order to solve the system.
  template <int dim>
  void ObstacleProblem<dim>::solve ()
  {
  std::cout << "  |> SOLVING THE SYSTEM (wait) ... " << std::endl;

    ReductionControl          reduction_control (1000, 1e-12, 1e-5);
    SolverCG<Vector<double>>  solver (reduction_control);

    precondition.initialize (bspline_system_matrix);

    solver.solve (bspline_system_matrix, bspline_solution, bspline_system_rhs, precondition);
    bspline_constraints.distribute (bspline_solution);

    cg_iter = reduction_control.last_step();    // [ ?? ]

    std::cout << "  |  Error: " << reduction_control.initial_value()
              << " -> " << reduction_control.last_value()
              << " in "
              <<  reduction_control.last_step()
              << " CG iterations."
              << std::endl;
  }

//-------------------------
/// Main problem class: This method construct and save the image output files. <br>
/*! In particular it prints on .vtk files the 3D plot of the function at every cycle step.
*/
  template <int dim>
  void ObstacleProblem<dim>::output_results (const unsigned int iteration)
  {
  std::cout << "  |- Saving/overwriting the .vtk result files." << std::endl;

    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);

    Vector<double> bspline_sol_dh(dof_handler.n_dofs());
    Vector<double> bspline_solution_vector(bspline_solution);

    iga_handler.transform_vector_into_fe_space(
      bspline_sol_dh,
      bspline_solution_vector);

    data_out.add_data_vector (bspline_sol_dh, "displacement");

    Vector<double> obstacle_vector(triangulation.n_active_cells());

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    unsigned int i=0;
    const Obstacle<dim> obstacle;
    for (; cell!=endc; ++cell, ++i)
      {
        obstacle_vector[i] = obstacle.value(cell->center());
      }

    data_out.add_data_vector (obstacle_vector, "obstacle");

    Vector<double> active_set_dh(dof_handler.n_dofs());
    iga_handler.transform_vector_into_fe_space(
      active_set_dh,
      active_set_vector);
    data_out.add_data_vector (active_set_dh, "active_set");

    Vector<double> bspline_contact_force_dh(dof_handler.n_dofs());
    Vector<double> bspline_contact_force_vector(bspline_contact_force);

    iga_handler.transform_vector_into_fe_space(
      bspline_contact_force_dh,
      bspline_contact_force_vector);

    data_out.add_data_vector (bspline_contact_force_dh, "lambda");

    data_out.build_patches ();

    //std::ofstream output_vtk ((std::string("output_") +
    //                           Utilities::int_to_string (iteration, 3) +
    //                           ".vtk").c_str ());

    std::string filename = "output_";
    filename += Utilities::int_to_string (iteration, 3);
    filename += ".vtk";
    std::string relpath = "RESULTS/" + filename;

    std::ofstream output_vtk (relpath);
    data_out.write_vtk (output_vtk);
  }

//-------------------------
/// Main problem class: This method process the solution in order to compute the L2 and the H1 norm.  <br>
/*! In addition it prepares the convergence table.
*/
  template <int dim>
  void ObstacleProblem<dim>::compute_error (unsigned int cycle)
  {
    std::cout << "\n - Computing the error..." << std::endl;
    Vector<float> difference_per_cell (triangulation.n_active_cells());

    Vector<double> bspline_sol_dh(dof_handler.n_dofs());
    Vector<double> bspline_solution_vector(bspline_solution);

    iga_handler.transform_vector_into_fe_space(
      bspline_sol_dh,
      bspline_solution_vector);

    VectorTools::integrate_difference (dof_handler, bspline_sol_dh,
                                       Solution<dim>(),
                                       difference_per_cell,
                                       QGauss<dim>(2*fe.degree+1),
                                       VectorTools::L2_norm);

    const double L2_error = difference_per_cell.l2_norm();

    VectorTools::integrate_difference (dof_handler, bspline_sol_dh,
                                       Solution<dim>(),
                                       difference_per_cell,
                                       QGauss<dim>(2*fe.degree+1),
                                       VectorTools::H1_seminorm);

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

//-------------------------
/// Main problem class: This is the function which has the top-level control over everything. 
  template <int dim>
  void ObstacleProblem<dim>::run (unsigned int cycle)
  {
	std::cout << " - Problem RUN" << std::endl;

    std::cout << " - Checking/creating the output directory \"RESULTS\". " << std::endl;
    std::string OutputFolderName = "RESULTS";
    std::filesystem::create_directories(OutputFolderName);

    make_grid();

    setup_system ();

    IndexSet active_set_old (active_set);

    std::cout << " > Starting the Newton iterations:\n" << std::endl;
    for (unsigned int iteration=0; iteration<=bspline_solution.size (); ++iteration)
      {
        std::cout << "  |[ Newton iteration: " << iteration << " ]" << std::endl;

        assemble_system ();

        if (iteration == 0)    // Initialization at iteration 0
          {
            bspline_complete_system_matrix.copy_from (bspline_system_matrix);
            bspline_complete_system_rhs = bspline_system_rhs;
          }

        solve ();
        update_solution_and_constraints ();
        output_results (iteration);

        if (active_set == active_set_old)
          {
            compute_error (cycle);    // [ IGA ONLY DIFFERENCE! ]
            break;
          }

        active_set_old = active_set;

        std::cout << std::endl;
      }
  }

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
      std::cout << "\n\n - Preparing the output converging table." << std::endl;
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
/**
* MAIN FUNCTION
* Entry point of the entire code.
*/
int main (int argc, char *argv[])
{
  std::cout << " " << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= THE CODE IS RUNNING ======================" << std::endl;
  std::cout << "============================================================" << std::endl;

  // DEFAULT VALUES:
  //--------------------------------------------------------
    unsigned int n_cycle       = 1;       // 5  (6)
    unsigned int degree        = 0;       // 0
    bool         h_refinement  = false;    // true
    bool         p_refinement  = false;   // false
    bool         k_refinement  = true;   // false
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
      using namespace dealii;
      using namespace Step41;
      deallog.depth_console (0);

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

		  // RUN: =========
		  std::cout << " - Assemble IgaHandler ...  " << std::endl;
		  Timings::Chrono Timer;		// NEW (START TIMER)
          Timer.start();				// NEW
            IgaHandler<2,2> iga_hd2(knots, mults, degree);
          Timer.stop();					// NEW (STOP TIMER)
          std::cout << "   (Time to assemble the IgaHandler: "
                    << Timer.wallTime()/1000000 << " seconds)" << std::endl;

          ObstacleProblem<2> obstacle_problem(iga_hd2, convergence_table);
          obstacle_problem.run (cycle);
          // END: =========

		  std::cout << " - CYCLE # " << cycle << " successfully ended!" << std::endl;
        }

    // Print the convergence table on screen and on a .txt file:
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

//====================================================
//====================================================
//====================================================


