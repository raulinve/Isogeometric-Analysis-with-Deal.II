/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2020 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */


// @sect3{Include files}

// As usual, the first few include files are already known, so we will not
// comment on them further.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/convergence_table.h>  // IGA new

#include <deal.II/lac/affine_constraints.h>  // ?
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>        // IGA
#include <deal.II/lac/sparse_direct.h>        // IGA
#include <deal.II/lac/precondition.h>         // IGA

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_bernstein.h>          // IGA
#include <deal.II/fe/mapping_fe_field.h>      // IGA

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>
#include <iostream>
#include <list>                               // IGA
#include <string>                             // NEW
#include <filesystem>                         // NEW

#include "grid_generator.h"                   // IGA HEADER
#include "iga_handler.h"                      // IGA HEADER


//namespace Step8
//{
  using namespace dealii;

  bool original_problem_8 = true;	// Default: "false" [SHOULD BE REMOVED NEXT]


//====================================================
  // @sect3{Right hand side values}

  template <int dim>
  void right_hand_side(const std::vector<Point<dim>> &points,
                       std::vector<Tensor<1, dim>> &  values)
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    Assert(dim >= 2, ExcNotImplemented());

    if (original_problem_8) {
		Point<dim> point_1, point_2;
		point_1(0) = 0.5;
		point_2(0) = -0.5;

		for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
		  {
		    // If <code>points[point_n]</code> is in a circle (sphere) of radius
		    // 0.2 around one of these points, then set the force in x-direction
		    // to one, otherwise to zero:
		    if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) ||
		        ((points[point_n] - point_2).norm_square() < 0.2 * 0.2))
		      values[point_n][0] = 1.0;
		    else
		      values[point_n][0] = 0.0;

		    // Likewise, if <code>points[point_n]</code> is in the vicinity of the
		    // origin, then set the y-force to one, otherwise to zero:
		    if (points[point_n].norm_square() < 0.2 * 0.2)
		      values[point_n][1] = 1.0;
		    else
		      values[point_n][1] = 0.0;
		  }
    }
    else
    {
        // [...]
    }
  }



//====================================================
//====================================================
  // @sect3{The <code>ElasticProblem</code> class template}

  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem(IgaHandler<dim,dim>  & iga_handler);
    ~ElasticProblem();
    void run(unsigned int cycle);

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle) const;

  IgaHandler<dim,dim>  &iga_handler;                    // IGA HEADER

  unsigned int         degree;                          // IGA
  unsigned int         cg_iter;                         // IGA  [ ?? used only in solve ]

    Triangulation<dim>  & triangulation;
    DoFHandler<dim>     & dof_handler;

    //FESystem<dim>   &fe;
    //FESystem<dim>(FE_Bernstein<dim>) &fe;
    //FESystem<dim>(FE_Bernstein<dim>(degree), dim)  &fe;   // XXX
    FE_Bernstein<dim>    &fe;   

    MappingFEField<dim>  *mappingfe;    

    AffineConstraints<double> bspline_constraints;

    SparsityPattern      sparsity_bspline;
    SparseMatrix<double> bspline_system_matrix;

    Vector<double> bspline_solution;
    Vector<double> bspline_system_rhs;
  };


  // @sect3{The <code>ElasticProblem</code> class implementation}

  // @sect4{ElasticProblem::ElasticProblem constructor}

  template <int dim>
  ElasticProblem<dim>::ElasticProblem(IgaHandler<dim,dim>  & iga_handler)
    :
    iga_handler(iga_handler),
    degree(iga_handler.degree),
    triangulation(iga_handler.tria),
    fe (iga_handler.fe),   // XXX .fe_sys
    dof_handler (iga_handler.dh),
    mappingfe(iga_handler.map_fe)
    //convergence_table(convergence_table)
  {
    std::cout << " - Initialization of the problem." << std::endl;
  }


  template <int dim>
  ElasticProblem<dim>::~ElasticProblem()
  {
    bspline_system_matrix.clear();
    if (mappingfe)
      delete mappingfe;
  }

  // @sect4{ElasticProblem::setup_system}


  template <int dim>
  void ElasticProblem<dim>::setup_system()
  {
    std::cout << " - Setup the system" << std::endl;

    dof_handler.distribute_dofs(fe);
    //solution.reinit(dof_handler.n_dofs());
    //system_rhs.reinit(dof_handler.n_dofs());
    bspline_solution.reinit (iga_handler.n_bspline);
    bspline_system_rhs.reinit (iga_handler.n_bspline);

    //----
    /*bspline_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, bspline_constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(dim),
                                             bspline_constraints);
    bspline_constraints.close();*/
    //----

    //DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    //DoFTools::make_sparsity_pattern(dof_handler,
    //                                dsp,
    //                                constraints,
    //                                /*keep_constrained_dofs = */ false);
    //sparsity_pattern.copy_from(dsp);
    //system_matrix.reinit(sparsity_pattern);

    DynamicSparsityPattern bspline_sp(iga_handler.n_bspline, iga_handler.n_bspline);
    iga_handler.make_sparsity_pattern (bspline_sp);
    sparsity_bspline.copy_from(bspline_sp);
    bspline_system_matrix.reinit(sparsity_bspline);


  }


  // @sect4{ElasticProblem::assemble_system}

  template <int dim>
  void ElasticProblem<dim>::assemble_system()
  {
  std::cout << " > ASSEMBLING THE SYSTEM (wait) ... " << std::endl;

    bspline_system_matrix = 0;
    bspline_system_rhs    = 0;

    //QGauss<dim> quadrature_formula(fe.degree + 1);
    const QGauss<dim>         quadrature_formula(fe.degree+1);
    //const RightHandSide<dim>  right_hand_side;

    FEValues<dim> fe_values(*mappingfe, fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // As was shown in previous examples as well, we need a place where to
    // store the values of the coefficients at all the quadrature points on a
    // cell. In the present situation, we have two coefficients, lambda and
    // mu.
    std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);

    // Well, we could as well have omitted the above two arrays since we will
    // use constant coefficients for both lambda and mu, which can be declared
    // like this. They both represent functions always returning the constant
    // value 1.0. Although we could omit the respective factors in the
    // assemblage of the matrix, we use them here for purpose of
    // demonstration.
    Functions::ConstantFunction<dim> lambda(1.), mu(1.);

    // Like the two constant functions above, we will call the function
    // right_hand_side just once per cell to make things simpler.
    std::vector<Tensor<1, dim>> rhs_values(n_q_points);

    // Now we can begin with the loop over all cells:
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    //for (const auto &cell : dof_handler.active_cell_iterators())
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs    = 0;

        fe_values.reinit(cell);

        // Next we get the values of the coefficients at the quadrature
        // points. Likewise for the right hand side:
        lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
        mu.value_list(fe_values.get_quadrature_points(), mu_values);
        right_hand_side(fe_values.get_quadrature_points(), rhs_values);


        // With this knowledge, we can assemble the local matrix
        // contributions:
        //for (const unsigned int i : fe_values.dof_indices())
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int component_i = fe.system_to_component_index(i).first;

            //for (const unsigned int j : fe_values.dof_indices())
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const unsigned int component_j = fe.system_to_component_index(j).first;

                //for (const unsigned int q_point : fe_values.quadrature_point_indices())
                for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                  {
                    cell_matrix(i, j) +=
                      (                                                  //
                        (fe_values.shape_grad(i, q_point)[component_i] * //
                         fe_values.shape_grad(j, q_point)[component_j] * //
                         lambda_values[q_point])                         //
                        +                                                //
                        (fe_values.shape_grad(i, q_point)[component_j] * //
                         fe_values.shape_grad(j, q_point)[component_i] * //
                         mu_values[q_point])                             //
                        +                                                //
                        ((component_i == component_j) ?        //
                           (fe_values.shape_grad(i, q_point) * //
                            fe_values.shape_grad(j, q_point) * //
                            mu_values[q_point]) :              //
                           0)                                  //
                        ) *                                    //
                      fe_values.JxW(q_point);                  //
                  }
              }
          }

        // Assembling the right hand side is also just as discussed in the
        // introduction:
        //for (const unsigned int i : fe_values.dof_indices())
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int component_i =
              fe.system_to_component_index(i).first;

            //for (const unsigned int q_point : fe_values.quadrature_point_indices())
            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              cell_rhs(i) += fe_values.shape_value(i, q_point) *
                             rhs_values[q_point][component_i] *
                             fe_values.JxW(q_point);
          }


        //cell->get_dof_indices(local_dof_indices);
        //constraints.distribute_local_to_global(
        //  cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);

        iga_handler.distribute_local_to_global(cell_matrix,
                                               cell_rhs,
                                               cell,
                                               bspline_system_matrix,
                                               bspline_system_rhs,
                                               bspline_constraints);

      }
  }



  // @sect4{ElasticProblem::solve}

  template <int dim>
  void ElasticProblem<dim>::solve()
  {
  std::cout << " > SOLVING THE SYSTEM (wait) ... " << std::endl;

    SolverControl            solver_control(100000, 1e-14);             // default: (1000, 1e-12)
    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(bspline_system_matrix, 1.2);

    cg.solve(bspline_system_matrix, bspline_solution, bspline_system_rhs, preconditioner);
    cg_iter = solver_control.last_step();    // [ ?? ]

    bspline_constraints.distribute(bspline_solution);

  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;
}


  // @sect4{ElasticProblem::refine_grid}
  template <int dim>
  void ElasticProblem<dim>::refine_grid()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1),
                                       {},
                                       bspline_solution,
                                       estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);

    triangulation.execute_coarsening_and_refinement();
  }


  // @sect4{ElasticProblem::output_results}

  template <int dim>
  void ElasticProblem<dim>::output_results(const unsigned int cycle) const
  {
  std::cout << " - Saving/overwriting the .vtk result files." << std::endl;

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    std::vector<std::string> solution_names;
    switch (dim)
      {
        case 1:
          solution_names.emplace_back("displacement");
          break;
        case 2:
          solution_names.emplace_back("x_displacement");
          solution_names.emplace_back("y_displacement");
          break;
        case 3:
          solution_names.emplace_back("x_displacement");
          solution_names.emplace_back("y_displacement");
          solution_names.emplace_back("z_displacement");
          break;
        default:
          Assert(false, ExcNotImplemented());
      }

    Vector<double> bspline_sol_dh(dof_handler.n_dofs());
    Vector<double> bspline_solution_vector(bspline_solution);
    iga_handler.transform_vector_into_fe_space(
      bspline_sol_dh,
      bspline_solution_vector);

    data_out.add_data_vector(bspline_sol_dh, solution_names);
    data_out.build_patches();

  std::string filename = "solution-";
  filename += ('0' + cycle);
  filename += ".vtk";
  std::string relpath = "RESULTS/" + filename;
  std::ofstream output(relpath);
  data_out.write_vtk(output);
  }



  // @sect4{ElasticProblem::run}

  template <int dim>
  void ElasticProblem<dim>::run(unsigned int cycle)
  {
	std::cout << " - Problem RUN (" << cycle << ")" << std::endl;

    std::cout << " - Checking/creating the output directory \"RESULTS\". " << std::endl;
    std::string OutputFolderName = "RESULTS";
    std::filesystem::create_directories(OutputFolderName);


    /*for (unsigned int cycle = 0; cycle < 8; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube(triangulation, -1, 1);
            triangulation.refine_global(4);
          }
        else
          refine_grid();

        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells() << std::endl;*/

        setup_system();

        std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl;

        assemble_system();
        solve();
        output_results(cycle);
      //}
  }
//} // namespace Step8







// @sect3{The <code>main</code> function}

int main(int argc, char *argv[])
{
  std::cout << " " << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= THE CODE IS RUNNING ======================" << std::endl;
  std::cout << "============================================================" << std::endl;

  // DEFAULT VALUES:
  //--------------------------------------------------------
    unsigned int n_cycle       = 3;       // 5  (6)
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
  //deallog.depth_console(0);
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  // Initializations before the loops:
  std::vector<unsigned int> subdivisions(1);
  //ConvergenceTable          convergence_table;

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
      ElasticProblem<2> elastic_problem_2d(iga_hd2);
      elastic_problem_2d.run(cycle);
      //=======================================

      // END: =========

	  std::cout << " - CYCLE # " << cycle << " successfully ended!" << std::endl;
    }

	//print_table(convergence_table, h_refinement, p_refinement, k_refinement);

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

  std::cout << " \n" << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= CODE ENDED CORRECTLY =====================" << std::endl;
  std::cout << "============================================================\n" << std::endl;

  return 0;
}
