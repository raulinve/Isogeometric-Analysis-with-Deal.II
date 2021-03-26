/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2011 - 2019 by the deal.II authors
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
 * Authors: Joerg Frohne, Texas A&M University and
 *                        University of Siegen, 2011, 2012
 *          Wolfgang Bangerth, Texas A&M University, 2012
 */
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
namespace Step41
{
  using namespace dealii;
  template <int dim>
  class ObstacleProblem
  {
  public:
    ObstacleProblem();
    void run();
  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void
         assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix);
    void update_solution_and_constraints();
    void solve();
    void output_results(const unsigned int iteration) const;
    Triangulation<dim>        triangulation;
    FE_Q<dim>                 fe;
    DoFHandler<dim>           dof_handler;
    AffineConstraints<double> constraints;
    IndexSet                  active_set;
    TrilinosWrappers::SparseMatrix system_matrix;
    TrilinosWrappers::SparseMatrix complete_system_matrix;
    TrilinosWrappers::MPI::Vector solution;
    TrilinosWrappers::MPI::Vector system_rhs;
    TrilinosWrappers::MPI::Vector complete_system_rhs;
    TrilinosWrappers::MPI::Vector diagonal_of_mass_matrix;
    TrilinosWrappers::MPI::Vector contact_force;
  };
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & /*p*/,
                         const unsigned int component = 0) const override
    {
      (void)component;
      AssertIndexRange(component, 1);
      return -10;
    }
  };
  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & /*p*/,
                         const unsigned int component = 0) const override
    {
      (void)component;
      AssertIndexRange(component, 1);
      return 0;
    }
  };
  template <int dim>
  class Obstacle : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override
    {
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
    }
  };
  template <int dim>
  ObstacleProblem<dim>::ObstacleProblem()
    : fe(1)
    , dof_handler(triangulation)
  {}
  template <int dim>
  void ObstacleProblem<dim>::make_grid()
  {
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(7);
    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "Total number of cells: " << triangulation.n_cells()
              << std::endl;
  }
  template <int dim>
  void ObstacleProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    active_set.set_size(dof_handler.n_dofs());
    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl
              << std::endl;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             BoundaryValues<dim>(),
                                             constraints);
    constraints.close();
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    system_matrix.reinit(dsp);
    complete_system_matrix.reinit(dsp);
    IndexSet solution_index_set = dof_handler.locally_owned_dofs();
    solution.reinit(solution_index_set, MPI_COMM_WORLD);
    system_rhs.reinit(solution_index_set, MPI_COMM_WORLD);
    complete_system_rhs.reinit(solution_index_set, MPI_COMM_WORLD);
    contact_force.reinit(solution_index_set, MPI_COMM_WORLD);
    TrilinosWrappers::SparseMatrix mass_matrix;
    mass_matrix.reinit(dsp);
    assemble_mass_matrix_diagonal(mass_matrix);
    diagonal_of_mass_matrix.reinit(solution_index_set);
    for (unsigned int j = 0; j < solution.size(); j++)
      diagonal_of_mass_matrix(j) = mass_matrix.diag_element(j);
  }
  template <int dim>
  void ObstacleProblem<dim>::assemble_system()
  {
    std::cout << "   Assembling system..." << std::endl;
    system_matrix = 0;
    system_rhs    = 0;
    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    RightHandSide<dim> right_hand_side;
    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs    = 0;
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) +=
                  (fe_values.shape_grad(i, q_point) *
                   fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));
              cell_rhs(i) +=
                (fe_values.shape_value(i, q_point) *
                 right_hand_side.value(fe_values.quadrature_point(q_point)) *
                 fe_values.JxW(q_point));
            }
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               cell_rhs,
                                               local_dof_indices,
                                               system_matrix,
                                               system_rhs,
                                               true);
      }
  }
  template <int dim>
  void ObstacleProblem<dim>::assemble_mass_matrix_diagonal(
    TrilinosWrappers::SparseMatrix &mass_matrix)
  {
    Assert(fe.degree == 1, ExcNotImplemented());
    const QTrapez<dim> quadrature_formula;
    FEValues<dim>      fe_values(fe,
                            quadrature_formula,
                            update_values | update_JxW_values);
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell_matrix = 0;
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_matrix(i, i) +=
              (fe_values.shape_value(i, q_point) *
               fe_values.shape_value(i, q_point) * fe_values.JxW(q_point));
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               local_dof_indices,
                                               mass_matrix);
      }
  }
  template <int dim>
  void ObstacleProblem<dim>::update_solution_and_constraints()
  {
    std::cout << "   Updating active set..." << std::endl;
    const double penalty_parameter = 100.0;
    TrilinosWrappers::MPI::Vector lambda(
      complete_index_set(dof_handler.n_dofs()));
    complete_system_matrix.residual(lambda, solution, complete_system_rhs);
    contact_force = lambda;
    contact_force.scale(diagonal_of_mass_matrix);
    contact_force *= -1;
    constraints.clear();
    active_set.clear();
    const Obstacle<dim> obstacle;
    std::vector<bool>   dof_touched(dof_handler.n_dofs(), false);
    for (const auto &cell : dof_handler.active_cell_iterators())
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          Assert(dof_handler.get_fe().dofs_per_cell ==
                   GeometryInfo<dim>::vertices_per_cell,
                 ExcNotImplemented());
          const unsigned int dof_index = cell->vertex_dof_index(v, 0);
          if (dof_touched[dof_index] == false)
            dof_touched[dof_index] = true;
          else
            continue;
          const double obstacle_value = obstacle.value(cell->vertex(v));
          const double solution_value = solution(dof_index);
          if (lambda(dof_index) + penalty_parameter *
                                    diagonal_of_mass_matrix(dof_index) *
                                    (solution_value - obstacle_value) <
              0)
            {
              active_set.add_index(dof_index);
              constraints.add_line(dof_index);
              constraints.set_inhomogeneity(dof_index, obstacle_value);
              solution(dof_index) = obstacle_value;
              lambda(dof_index) = 0;
            }
        }
    std::cout << "      Size of active set: " << active_set.n_elements()
              << std::endl;
    std::cout << "   Residual of the non-contact part of the system: "
              << lambda.l2_norm() << std::endl;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             BoundaryValues<dim>(),
                                             constraints);
    constraints.close();
  }
  template <int dim>
  void ObstacleProblem<dim>::solve()
  {
    std::cout << "   Solving system..." << std::endl;
    ReductionControl                        reduction_control(100, 1e-12, 1e-3);
    SolverCG<TrilinosWrappers::MPI::Vector> solver(reduction_control);
    TrilinosWrappers::PreconditionAMG       precondition;
    precondition.initialize(system_matrix);
    solver.solve(system_matrix, solution, system_rhs, precondition);
    constraints.distribute(solution);
    std::cout << "      Error: " << reduction_control.initial_value() << " -> "
              << reduction_control.last_value() << " in "
              << reduction_control.last_step() << " CG iterations."
              << std::endl;
  }
  template <int dim>
  void ObstacleProblem<dim>::output_results(const unsigned int iteration) const
  {
    std::cout << "   Writing graphical output..." << std::endl;
    TrilinosWrappers::MPI::Vector active_set_vector(
      dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
    for (const auto index : active_set)
      active_set_vector[index] = 1.;
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "displacement");
    data_out.add_data_vector(active_set_vector, "active_set");
    data_out.add_data_vector(contact_force, "lambda");
    data_out.build_patches();
    std::ofstream output_vtk("output_" +
                             Utilities::int_to_string(iteration, 3) + ".vtk");
    data_out.write_vtk(output_vtk);
  }
  template <int dim>
  void ObstacleProblem<dim>::run()
  {
    make_grid();
    setup_system();
    IndexSet active_set_old(active_set);
    for (unsigned int iteration = 0; iteration <= solution.size(); ++iteration)
      {
        std::cout << "Newton iteration " << iteration << std::endl;
        assemble_system();
        if (iteration == 0)
          {
            complete_system_matrix.copy_from(system_matrix);
            complete_system_rhs = system_rhs;
          }
        solve();
        update_solution_and_constraints();
        output_results(iteration);
        if (active_set == active_set_old)
          break;
        active_set_old = active_set;
        std::cout << std::endl;
      }
  }
} // namespace Step41
int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step41;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
      AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
                  ExcMessage(
                    "This program can only be run in serial, use ./step-41"));
      ObstacleProblem<2> obstacle_problem;
      obstacle_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
