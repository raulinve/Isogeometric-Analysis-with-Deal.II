/**
* @file         beam.cc
* @brief        beam code
* @detail       Isogeometric Analysis for the solution of an elastic equation using the deal.II library.
* @author       Raul Invernizzi.
*/
// @include      string.h

/*! \mainpage MAIN PAGE
*
* Isogeometric Analysis classes for the deal.II library.
* 
* <hr>
* This code is used to introduce the Isogeometric Analysis (IGA) framework 
* into the standard deal.II Finite Element (FE) library.
* 
* In this specific example named "beam" a bezier polynomia approximation of 
* the elastic equation is computed on a structure with a Cook's beam geoemtry.  <br> 
* 
*
* Please refere to the most recent version available on GitHub at the following link:  
* <a href="https://github.com/raulinve/Isogeometric-Analysis-with-Deal.II">
* Isogeometric-Analysis-with-Deal.II</a>. \n\n
*
* <hr>
* <b>Usage:</b>  
*
* The executable takes 0 or 5 arguments to run, respectively: <br> 
* <pre>
* <tt><b>    ./beam  </b></tt>
*     or
* <tt><b>    ./beam  FE_TYPE  QUADRATURE  DEG  CYCLE_START  CYCLE_STOP </b></tt>
* </pre>
* 
* The accepted arguments of the program are: <br> 
* <pre>
*    FE_TYPE:
*        bernstein [default]
*        lagrange
*        lobatto
*    QUADRATURE:
*        legendre [default]
*        lobatto
*    DEG:         (es: 1) degree of the finite element space
*    CYCLE_START: (es: 0) initial refinement of the grid
*    CYCLE_STOP:  (es: 5) final refinement of the grid
* </pre>
* 
* For example, the default command is: <tt><b>./beam bernstein legendre 1 0 5 </b></tt>. \n
* 
* <hr>
* <b>Classes & Methods:</b>  
*
* int <b>main</b> (int argc, char **argv)<br> 
* 
* class <b>"BoundaryValues"</b>:<br> 
* class BoundaryValues : public Function<dim>(2) {}
* double BoundaryValues<dim>::value  
* BoundaryValues::value (solution at the boundary) :  
* 
* class <b>"RHSValues"</b>:<br> 
* class  RHSValues : public Function<dim>(2) {}
* double  RHSValues<dim>::value  
* 
* class <b>"ElasticProblem"</b>:<br> 
* class ElasticProblem 
* ElasticProblem<dim>::ElasticProblem  
* ElasticProblem<dim>::~ElasticProblem ()  
* void ElasticProblem<dim>::make_grid ()  
* void ElasticProblem<dim>::setup_system ()  
* void ElasticProblem<dim>::assemble_system ()  
* void ElasticProblem<dim>::solve ()  
* void ElasticProblem<dim>::output_results  
* void ElasticProblem<dim>::run ()  
* 
* <hr>
*/


/*! ---------------------------------------------------------------------
* Copyright (C) 1999 - 2020 by the deal.II authors
* Copyright (C) 2021 by Raul Invernizzi
*
* This file has been modified from the example program step-8 of the
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



// Include files:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_ilu.h>     // IGA

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>    // [ ADD ]
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>       // vector-valued finite elements
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe.h>              // IGA
#include <deal.II/fe/fe_nothing.h>      // IGA
#include <deal.II/fe/fe_bernstein.h>    // IGA
#include <deal.II/fe/mapping_q1.h>      // ADD NEW
#include <deal.II/fe/mapping.h>         // ADD NEW

#include <fstream>
#include <iostream>
#include <filesystem>                   // IGA


// Variables used in the development to activate 
// sections of code, set as "false".
bool IMPOSE_FORCE_RHS = false;
bool CHECK_TIP_DISP   = false;

bool PRINT_RESULTS = false;  // <- BETTER TO USE ONLY FOR DEBUGGING


using namespace dealii;
namespace Beam
{

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
			BoundaryValues () : Function<dim>(2) {}

			virtual double value (const Point<dim>   &p,
					                  const unsigned int  component = 0) const 
			{
				return 0*p(0)*component;    // <- [ USELESS - Just to avoid warnings of unused params ]
			}
	};


//====================================================
	/**
	* @class   right_hand_side
	* @brief   Class used to define the external forces applied on the structure.
	*
	* This class implement the forcs that are applied to the structure. 
  * These are stored in the so-called Righ Hand Side (RHS) vector.
	*/

  template <int dim>
  void right_hand_side(const std::vector<Point<dim>> & points,
                       std::vector<Tensor<1, dim>>   & values)
  { 
    /*
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    Assert(dim >= 2, ExcNotImplemented());

		// Points of application:
		Point<dim> point_1;
		point_1(0) =  46;   point_1(1) =  50;

		double radius = 2;     // Radius of the force footprint
		double force  = 0;     // Modulus of the forces [-0.1]

		for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
			{ 
			  if (((points[point_n] - point_1).norm_square() < radius * radius))
			    values[point_n][1] = force;
			  else
			    values[point_n][1] = 0.0;
			}
    */
    // DELETE THE FOLLOWING LINE:
    if(0) {values[0][0]=0*points[0](0);}  // <- [ USELESS - Just to avoid warnings of unused params ]
    
  }



	/**
	* @class   RHSValues
	* @brief   Class used to define the external forces applied on the structure.
	*
	* This class implement the forcs that are applied to the structure. 
  * It is an alternative with respect to the class right_hand_side,
  * The definition of the external forces in the case of Bernstein 
  * polynomia must be done using this class.
  * These forces are stored in the so-called Righ Hand Side (RHS) vector.
	*/
	template <int dim>
	class RHSValues : public Function<dim>
	{
		public:
			RHSValues () : Function<dim>(2) {}

			virtual double value (const Point<dim>   &p,
				                    const unsigned int  component = 0) const 
			{
				if (component == 1)                    // If Y component
				  return 1;                            // Apply the force
				else
				  return 0*p(0)*component; // <- [ USELESS - Just to avoid warnings of unused params ]
			}
	};

  double force_intensity(const unsigned int cycle,
                         const unsigned int degree) {
        double total_resultant_force      = 1.;      //[N]
        double right_side_size            = 16.;     //[mm]

				double global_force_on_right_side = total_resultant_force/right_side_size;  //[N/mm]
        double right_side_num_of_cells    = pow(2.,cycle);
  			double right_side_cell_size       = right_side_size / right_side_num_of_cells;

				double local_force_on_right_side = global_force_on_right_side*(right_side_cell_size)/2./degree;

				// Rekap:
        std::cout << "     TOTAL FORCE APPLIED           : " << global_force_on_right_side*right_side_size << " N" << std::endl;
				std::cout << "     Force magnitude (global)      : " << global_force_on_right_side << " N/mm" << std::endl;
			  std::cout << "     Number of cells on right side : " << right_side_num_of_cells    << " " << std::endl;
			  std::cout << "     Cell dimension                : " << right_side_cell_size       << " mm" << std::endl;
				std::cout << "     Force magnitude (local nodal) : " << local_force_on_right_side  << " N" << std::endl;

			  return local_force_on_right_side;
	}

//====================================================
/**
* @class   ElasticProblem
* @brief   Main problem class
* 
* This class implement and collect every aspect of the poisson problem.
*
* @note This code admits some parameters from the
* command line, and they are passed to the constructor of the
* problem as additional arguments.
*
* @note The FESystem class is used in order to construct a vector-valued
* finite element that is composed of several scalar finite elements.
* we pass it the finite element of which we would like to compose the
* system of, and how often it shall be repeated
*
* The accepted arguments of the program are:
* fe_name          : bernstein, lagrange, lobatto
* quadrature_name  : legendre, lobatto
* degree           : degree of the finite element space
* n_cycles_up  : initial refinement of the grid
* n_cycles_down  : final refinement of the grid
* 
*/
  template <int dim>
  class ElasticProblem
  {
  public:
    ElasticProblem(const std::string  fe_name,           // IGA: Addeed parameters
                   const std::string  quadrature_name,
                   const unsigned int degree,
                   const unsigned int n_cycles_low,
                   const unsigned int n_cycles_up);
    ~ElasticProblem();                                   // IGA

    void run();

  private:
    void make_grid(const unsigned int ref_cycle);
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
    FESystem<dim>        *fe;                // IGA: from fe(1) to *fe
    DoFHandler<dim>      dof_handler;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    std::map<types::global_dof_index, double> boundary_values;
    unsigned int tip_low, tip_mid, tip_top;

    Quadrature<dim>      matrix_quad;        // IGA
    Quadrature<dim>      error_quad;         // IGA
    Quadrature<dim-1>    boundary_quad;      // IGA

    // Plane Stress Formulation:
    double E_mod        = 2.; //30.0;     // Elastic modulus (1 MPa = 1 N/mm²)
    double poisson      = 1.0/3.0;  // Poisson coefficient (~ 0.0 - 0.5)

    double global_scale = 1.e-3;
  };


/// Main problem class: Constructor
  template <int dim>
  ElasticProblem<dim>::ElasticProblem(const std::string  fe_name,
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
    fe              (NULL),                 //fe(FE_Q<dim>(1), dim)
    dof_handler     (triangulation)
  {
    // IGA: This section is entirely new.
    if (quadrature_name == "legendre")      // [DEFAULT]
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
      fe = new FESystem<dim>(FE_Bernstein<dim>(degree), dim);
    else if (fe_name == "lagrange")         // [DEFAULT]
      fe = new FESystem<dim>(FE_Q<dim>(degree), dim);
    else if (fe_name == "lobatto")
      fe = new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(degree+1)), dim);
    else
	  AssertThrow(false, ExcMessage("FE not supported"));

  }


/// Main problem class: Destructor
	template <int dim> 
	ElasticProblem<dim>::~ElasticProblem ()                  // IGA ONLY
	{
	  dof_handler.clear();
	  if (fe)
		delete fe;
	}


/// Additional function: Transformation function to get the wanted beam-shape from the standard rectangle.
  template <int dim>
  Point<dim> grid_y_transform (const Point<dim> &pt_in)    // [ ADD ]
  {
    const double &x = pt_in[0];
    const double &y = pt_in[1];

    const double y_upper = 44.0 + (16.0/48.0)*x; // Line defining upper edge of beam
    const double y_lower =  0.0 + (44.0/48.0)*x; // Line defining lower edge of beam
    const double theta = y/44.0; // Fraction of height along left side of beam
    const double y_transform = (1-theta)*y_lower + theta*y_upper; // Final transformation

    Point<dim> pt_out = pt_in;
    pt_out[1] = y_transform;

    return pt_out;
  }


/// Main problem class: Method used to generate the geometry and the grid.
/*! NOTE: In this method, the different edges are also named with specific tags. 
*   In particular:
*     [1]  = Left side (clamped), 
*     [2]  = Upper and Lower edges (free), 
*     [11] = Right side (free, with forces)
*/
  template <int dim>
  void ElasticProblem<dim>::make_grid(const unsigned int ref_cycle)
  {
		// Here we create the triangulation of the domain, for which we choose 
		// a scaled an anisotripically discretised rectangle which is subsequently
		// transformed into the correct of the Cook cantilever.
		// Each relevant boundary face is then given a boundary ID number.

		// Divide the beam, but only along the x- and y-coordinate directions
		unsigned int elements_per_edge = pow(2,ref_cycle);
		std::vector< unsigned int > repetitions(dim, elements_per_edge);
		// Only allow one element through the thickness
		// (modelling a plane strain condition)
		if (dim == 3)
			repetitions[dim-1] = 1;

		const Point<dim> bottom_left = (dim == 3 ? Point<dim>( 0.0,  0.0, -0.5) : Point<dim>( 0.0,  0.0));
		const Point<dim> top_right   = (dim == 3 ? Point<dim>(48.0, 44.0,  0.5) : Point<dim>(48.0, 44.0));

		GridGenerator::subdivided_hyper_rectangle(triangulation,
			                                        repetitions,
			                                        bottom_left,
			                                        top_right);

    // CUSTOM BOUNDARY ID:
		// Since we wish to apply a Neumann BC to the right-hand surface, we
		// must find the cell faces in this part of the domain and mark them with
		// a distinct boundary ID number.  The faces we are looking for are on the
		// +x surface and will get boundary ID 11.
		// Dirichlet boundaries exist on the left-hand face of the beam (this fixed
		// boundary will get ID 1) and on the +Z and -Z faces (which correspond to
		// ID 2 and we will use to impose the plane strain condition)
		const double tol_boundary = 1e-6;
		typename Triangulation<dim>::active_cell_iterator cell =
			       triangulation.begin_active(), endc = triangulation.end();
		for (; cell != endc; ++cell)
			for (unsigned int face = 0;
			     face < GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary() == true)
			    {
			      if (std::abs(cell->face(face)->center()[0] - 0.0) < tol_boundary)
			        cell->face(face)->set_boundary_id(1);       // -X faces
			      else if (std::abs(cell->face(face)->center()[0] - 48.0) < tol_boundary)
			        cell->face(face)->set_boundary_id(11);      // +X faces
			      else if (std::abs(std::abs(cell->face(face)->center()[0]) - 0.5) < tol_boundary)
			        cell->face(face)->set_boundary_id(2);       // +Z and -Z faces
			    }

    // TRANSFORMATION OF THE GEOMETRY:   [STD SQUARE] -> [BEAM SHAPE]
		// Transform the hyper-rectangle into the beam shape
		GridTools::transform(&grid_y_transform<dim>, triangulation);

    // Optinal: Scaling
		// GridTools::scale(1e-3, triangulation);             // parameters.scale
	  
    // Optional: Volume
		// vol_reference = GridTools::volume(triangulation);
		// vol_current = vol_reference;
		// std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;
  }


/// Main problem class: Method used to setup the matrix system.
  template <int dim>
  void ElasticProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(*fe);                            // IGA: "*"
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

	  // IMPOSE CONSTRAINTS: -------------
    constraints.clear();

    const int boundary_id = 1;   // SELECT BOUNDARY TO FIX
    // [1] = Left side clamped, [2] = Upper and Lower edges, [11] = Right side free

    // NOTE: "interpolate" method works fine with std Lagrange polynomial, BUT FAILS with Bernstein. 
    /*const FEValuesExtractors::Vector u_fe(0);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_id,
                                             ZeroFunction<dim>(dim),
                                             constraints);
                                             fe->component_mask(u_fe));  */

    // NOTE: This is why it is needed to replece it with the following "project" method.
    std::map< types::boundary_id, const Function< dim > *>  dirichlet_boundary;        // IGA
    BoundaryValues<dim> boundary_funct;                                                // IGA
    dirichlet_boundary[boundary_id] = &boundary_funct;                                 // IGA

    VectorTools::project_boundary_values (dof_handler, 
                                          dirichlet_boundary,    // ZeroFunction<dim>(dim),
                                          boundary_quad,
                                          constraints);
    constraints.close();
	  // ---------------------------------

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);

    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);


    // IMPOSE RHS FORCES: --------------
    RHSValues<dim> rhs_function;
    std::map<types::boundary_id, const Function<dim, double> *> bc_map;
    bc_map[11] = &rhs_function;   // [11] referes to the free right edge
    boundary_values.clear();
    VectorTools::project_boundary_values(dof_handler,
                                         bc_map,
                                         boundary_quad,
                                         boundary_values,
                                         {0,1}); // component y.
    // NOTE: With Bernstein "interpolate does not work. Must use project"
    // VectorTools::interpolate_boundary_values(dof_handler,
    //                                      11,
    //                                      rhs_function,
    //                                      boundary_values,
    //                                      {false, true});
	  // ---------------------------------
  }


/// Main problem class: Method used to assemble the main system.
  template <int dim>
  void ElasticProblem<dim>::assemble_system()
  {
    std::cout << "   > ASSEMBLING THE SYSTEM (wait) ... " << std::endl;
    // QGauss<dim> quadrature_formula(fe.degree + 1);           // IGA deleted

    FEValues<dim> fe_values(*fe, matrix_quad,     // IGA: "*", "quadrature_formula" to "matrix_quad"
                            update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe->dofs_per_cell;       // IGA: "->"
    const unsigned int   n_q_points    = matrix_quad.size();      // IGA from "quadrature_formula.size();"

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs    (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Addition in case variable coefficients are needed:
    std::vector<double> lambda_values (n_q_points);
    std::vector<double> mu_values     (n_q_points); 

    //[N/mm], (1 MPa = 1 N/mm²)
    Functions::ConstantFunction<dim> lambda(E_mod*poisson/((1.+poisson)*(1.0-2.0*poisson)));   // 1st Lame Const (52.5)
    Functions::ConstantFunction<dim> mu(E_mod/(2.0*(1.0+poisson)));      // Shear modulus (=G) (26.25)

    // Like the two constant functions above, we will call the function
    // right_hand_side just once per cell to make things simpler.
    std::vector<Tensor<1, dim>>  rhs_values(n_q_points);


    typename DoFHandler<dim>::active_cell_iterator              // IGA new
    cell = dof_handler.begin_active(),                          // IGA new
    endc = dof_handler.end();                                   // IGA new
    // Now we can begin with the loop over all cells:
    for (; cell!=endc; ++cell)                                  // IGA mod
    // for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs    = 0;

        // Get the values of the coefficients at the quadrature points.
        lambda. value_list(fe_values.get_quadrature_points(), lambda_values); 
        mu.     value_list(fe_values.get_quadrature_points(), mu_values);
        right_hand_side   (fe_values.get_quadrature_points(), rhs_values);   // XXX

        // Then assemble the entries of the local stiffness matrix and right
        // hand side vector.
        // We can assemble the local matrix contributions:
        for (unsigned int i=0; i<dofs_per_cell; ++i)            // IGA mod
        // for (const unsigned int i : fe_values.dof_indices())
          {
            const unsigned int component_i = fe->system_to_component_index(i).first;     // IGA: "->"

            for (unsigned int j=0; j<dofs_per_cell; ++j)        // IGA mod
            // for (const unsigned int j : fe_values.dof_indices())
              {
                const unsigned int component_j =fe->system_to_component_index(j).first;  // IGA: "->"

                for (unsigned int q_point=0; q_point<n_q_points; ++q_point)    // IGA mod
                // for (const unsigned int q_point : fe_values.quadrature_point_indices()))
                  {
                    cell_matrix(i, j) +=
                      // CHECK DEAL.II "STEP-8" FOR DETAILS:
                      (                                                  //
                        (fe_values.shape_grad(i, q_point)[component_i] * //
                         fe_values.shape_grad(j, q_point)[component_j] * //
                         lambda_values[q_point])                         //
                        +                                                //
                        (fe_values.shape_grad(i, q_point)[component_j] * //
                         fe_values.shape_grad(j, q_point)[component_i] * //
                         mu_values[q_point])                             //
                        +                                                //
                        // CHECK DEAL.II "STEP-8" FOR DETAILS:
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

        // Assembling the right hand side:
        // STANDARD mode: Use the definition on the quadrature points:
        // NOTE: De-comment for usiing with the class "right_hand_side"
        // -----------------------
        //for (unsigned int i=0; i<dofs_per_cell; ++i)            // IGA mod
        //  {
        //    const unsigned int component_i = fe->system_to_component_index(i).first;     // IGA: "->"
        //    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)    // IGA mod
        //      cell_rhs(i) += fe_values.shape_value(i, q_point) *
        //                     rhs_values[q_point][component_i] *
        //                     fe_values.JxW(q_point);
        //  }
        // -----------------------

        // The transfer from local degrees of freedom into the global matrix
        // and right hand side vector does not depend on the equation under
        // consideration, and is thus the same as in all previous examples.
        cell -> get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               cell_rhs, 
                                               local_dof_indices, 
                                               system_matrix, 
                                               system_rhs);
      }

      // IMPOSE RHS FORCES: --------------
      //std::cout << "Print the System RHS:" << std::endl;
      //std::cout << " size 1: " << boundary_values.size() << "\n"; 
      unsigned int counter = 0;
      for (auto x : boundary_values)
        {
          //std::cout << x.first << ", " << x.second << std::endl;
          if (x.second != 0) {   // Insert the forces only if different from zero:
            counter++;
            if (counter == 1) {tip_low = x.first;}
            if (counter == boundary_values.size()/2) {tip_top = x.first;}
            if (counter == (boundary_values.size()/2-1)/2+1) {tip_mid = x.first;}
            system_rhs[x.first] = x.second;
            if (counter>1 && counter<boundary_values.size()/2) {
            	system_rhs[x.first] = 2.*x.second;
            }
          }
        }
	    // ---------------------------------

  }


/// Main problem class: This method is called in order to solve the system using the CG solver.
  template <int dim>
  void ElasticProblem<dim>::solve()
  {
  std::cout << "   > SOLVING THE SYSTEM (wait) ... " << std::endl;

    SolverControl             solver_control(100000, 1e-14);   // default: (1000, 1e-12)
    SolverCG<Vector<double>>  cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    std::cout << "   Memory consumption " << system_matrix.memory_consumption()
              << " bytes" << std::endl;

    cg.solve(system_matrix, 
             solution, 
             system_rhs, 
             preconditioner);

    constraints.distribute(solution);

    std::cout << "   " << solver_control.last_step()
              << " CG iterations needed to obtain convergence." << std::endl;

    cg_iter = solver_control.last_step();
  }


/// Main problem class: This method construct and save the image output files. <br>
/*! In particular it prints on .vtk files the 3D plot of the function 
*   at every cycle step.
*   The only difference is that the solution function is vector
*   valued. The DataOut class takes care of this automatically, but we have
*   to give each component of the solution vector a different name.
*   
*   To do this, the DataOut::add_vector() function wants a vector of
*   strings. Since the number of components is the same as the number
*   of dimensions we are working in, we use the <code>switch</code>
*   statement below.
*/
  template <int dim>
  void ElasticProblem<dim>::output_results(const unsigned int cycle) const
  {
    std::cout << "   Saving/overwriting the .vtk result files." << std::endl;
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    // After setting up the names for the different components of the
    // solution vector, we can add the solution vector to the list of
    // data vectors scheduled for output.
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
    data_out.add_data_vector(solution, solution_names);
    
    data_out.build_patches();

	std::string filename = "solution-" + std::to_string(dim) + "d-";
	filename += ('0' + cycle);
	filename += ".vtk";
    std::string relpath = "RESULTS/" + filename;    // ADD
    std::ofstream output (relpath);
    data_out.write_vtk(output);

  }


/// Main problem class: This is the function which has the top-level control over everything. 
  template <int dim>
  void ElasticProblem<dim>::run()
  {
    // Initial rekap (on-screen print) -------
    std::cout << " \n                    polyn     quad  deg refinements" << std::endl;
    std::cout << "SELECTED OPTIONS : " << fe_name << " "
                                       << quadrature_name << " "
                                       << degree << " "
                                       << n_cycles_low << " "
                                       << n_cycles_up << "\n" << std::endl;
                                       
    std::cout << "Checking/creating the output directory \"RESULTS\". " << std::endl;
    std::string OutputFolderName = "RESULTS";
    std::filesystem::create_directories(OutputFolderName);

    std::cout << "Solving problem in " << dim << " space dimensions. \n\n" << std::endl;   


    // MAIN CYCLE OVER THE CYCLES:
    // ===================================
    for (unsigned int cycle = n_cycles_low; cycle < n_cycles_up+1; ++cycle)  // IGA cycles
      {
        std::cout << "\n\nCycle " << cycle << " ( " << n_cycles_low << " to " << n_cycles_up << " ):" << std::endl;
        
        if (cycle == n_cycles_low)
          {
            make_grid(n_cycles_low);          // IGA ADD (+custom geometry)
          }
        else 
          {
            triangulation.refine_global(1);   // IGA: constant refinement
          }

		    std::cout << "   Grid of elements:      "
				          << sqrt(triangulation.n_active_cells()) << "x" << sqrt(triangulation.n_active_cells()) << std::endl
				          << "   Total number of cells: " << triangulation.n_cells() << std::endl;

		    setup_system();         // [MAIN METHOD]

		    std::cout << "   Total number of dofs:  " << dof_handler.n_dofs() << std::endl;

		    assemble_system();      // [MAIN METHOD]


				// The following part, if activeted, overwrite the rhs vector.
				// NOTE: It does not work with bernstein polynomia!
				// ---------------------------------
				if (IMPOSE_FORCE_RHS) 
				{
				  /*
					std::cout << "   - Applying rhs forces" << std::endl;

					// 01. Application points:
					// The beam tip is at x=48, y=44 up to y=44+16=60.
					std::vector<Point<dim>> application_pts;
					double refin_incr = 1./pow(2.,cycle);
					for (double y = 44.0; y<=44.0+16.0; y+=refin_incr*16.0) {
						application_pts.push_back(Point<dim>({48.0, y}));
				    if(y!=44.0 || y!=44.0+16.0) { 
							application_pts.push_back(Point<dim>({48.0, y}));
				    }
					}
					double tol = 1.0/1000000;
					MappingQ1< dim > stMap;

					// 02. Define force intensity:
					double force_mod_y = 6.25;//00;  //global_scale;  //[N/mm], (1 MPa = 1 N/mm²)
					double force_app_y = force_mod_y/2.*(refin_incr*16.);  //application_pts.size()*

					// Rekap:
				  std::cout << "     Cell dimension:  " << refin_incr*16 << std::endl;
					std::cout << "     Force magnitude: " << force_app_y << " ( = " << force_mod_y << "/" << application_pts.size() << " )" << std::endl;
					//std::cout << "     Force application nodes (x y):" << std::endl;  
					//std::cout << "      ";
					//int cnt = 0;
					//for(auto i : application_pts) {std::cout << "(" << i << ")  ";  ++cnt;
					//                               if(cnt==5){std::cout << std::endl; std::cout << "      "; cnt=0;} }  std::cout << std::endl;

					// FORCES in X-direction:
					std::map<types::global_dof_index, Point<dim>> support_points_x;
					std::vector< bool > x_c{true, false};
					ComponentMask cmask_x(x_c);

				  // NOTE: THE FOLLOWING FUNCTION DOES NOT WORK WITH BERNSTEIN!
					DoFTools::map_dofs_to_support_points(stMap,
										           			 dof_handler,
										           			 support_points_x,
										           			 cmask_x);

					for (const auto p : support_points_x)
					{
						//for (auto pt : application_pts) {
						//   if (p.second.distance(pt) < tol) { system_rhs[p.first] = 0; }
						//}
						system_rhs[p.first] = 0;
					}

					// FORCES in Y-direction:
					std::map<types::global_dof_index, Point<dim>> support_points_y;
					std::vector< bool > y_c{false, true};
					ComponentMask cmask_y(y_c);
					DoFTools::map_dofs_to_support_points(stMap,
										          			 dof_handler,
										          			 support_points_y,
										          			 cmask_y);
					for (const auto p : support_points_y)
					{
						for (auto pt : application_pts) {
								//if (support_points_y[i].distance(pt) < tol) { system_rhs[i] = force_app_y; }
								if (p.second.distance(pt) < tol) { system_rhs[p.first] = force_app_y; }
						}
					}
				  */
				}
				// ---------------------------------
        // Print the rhs vector:
        //if (PRINT_RHS) {
		    //  std::cout << "  Print full rhs (IDs only!):" << std::endl;
		    //  unsigned int node_num = 0;
				//	for (auto j : system_rhs) {
		    //    std::cout << " " << node_num++ << "\t: " << j << std::endl;
				//	}
        //}

        // Compute the real force per node and apply to the rhs:
				double force_int = force_intensity(cycle, degree);
			  for (unsigned int i=0; i<system_rhs.size(); ++i) { system_rhs[i] *= force_int; }


				solve();                // [MAIN METHOD]

	      std::cout << "\n   Print the RHS and the SOLUTION at the tip:" << std::endl;
				std::cout << "          \t: " << "ID" << "\t: " << "RHS" << "\t: " << "SOL" << std::endl;
				std::cout << "   tip_top\t: " << tip_top << "\t: " << system_rhs[tip_top] << "\t: " << solution[tip_top] << std::endl;
				std::cout << "   tip_mid\t: " << tip_mid << "\t: " << system_rhs[tip_mid] << "\t: " << solution[tip_mid] << std::endl;
				std::cout << "   tip_low\t: " << tip_low << "\t: " << system_rhs[tip_low] << "\t: " << solution[tip_low] << std::endl;
				std::cout << std::endl;

        if (PRINT_RESULTS) {
		      std::cout << "\n   Print the full RHS and the SOLUTION vector:" << std::endl;
					std::cout << "   " << "ID" << "\t: " << "RHS" << "\t: " << "SOL" << std::endl;
		      unsigned int node_num = 0;
					for (unsigned int i = 0; i<solution.size(); i++) {
		        std::cout << "   " << node_num++ << "\t: " << system_rhs[i] << "\t: " << solution[i] << std::endl;
					}
        }

				// ---------------------------------
				// The following part, if activeted, prints the tip vertical displacement on screen.
				// NOTE: It does not work with bernstein polynomia!
				if (CHECK_TIP_DISP) 
				{ 
				  /*
					// Check tip displacement:
					std::cout << "   - Check tip displacement" << std::endl;
					std::cout << "     Tip vertical displacement:" << std::endl;
					// The beam tip is at x=48, y=44 up to y=44+16=60.
					std::map<types::global_dof_index, Point<dim>> sup_pts;
					//std::vector< bool > y_c{false, true};
					//ComponentMask cmask_y(y_c);
						//++++++
						 double tol = 1.0/1000000;
						 MappingQ1< dim > stMap;
						 std::map<types::global_dof_index, Point<dim>> support_points_y;
						 std::vector< bool > y_c{false, true};
						 ComponentMask cmask_y(y_c);
						//++++++
					DoFTools::map_dofs_to_support_points(stMap, dof_handler, sup_pts, cmask_y);
					for (const auto p : sup_pts)
					{
						if (p.second.distance(Point<dim>({48, 44})) < tol) { std::cout << "     (low tip): " << solution[p.first]; }
						if (p.second.distance(Point<dim>({48, 52})) < tol) { std::cout << "  (mid tip): " << solution[p.first]; }
						if (p.second.distance(Point<dim>({48, 60})) < tol) { std::cout << "  (top tip): " << solution[p.first] << std::endl; }
					}
				  */
				}
				// ---------------------------------

        output_results(cycle);      // [MAIN METHOD]
    }
    // ===================================
  }




} // namespace Beam




//====================================================
/**
* MAIN FUNCTION
* Entry point of the entire code.
*/
int main(int argc, char **argv)
{
  std::cout << " " << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= THE CODE IS RUNNING ======================" << std::endl;
  std::cout << "============================================================" << std::endl;

  // DEFAULT VALUES:
  //--------------------------------------------------------
  char fe_name[]             = "lagrange"; // “lagrange”
  char quad_name[]           = "legendre";  // “legendre”
  unsigned int degree        = 2;           // 1
  unsigned int n_cycles_down = 0;           // 0   (Pieces: 0>1, 1>2, 2>4, 3>8, 4>16, 5>32, 6>64)
  unsigned int n_cycles_up   = 4;           // 2
  //--------------------------------------------------------
  // EXAMPLES:
  // ./beam lagrange legendre 1 0 4
  // ./beam bernstein legendre 1 0 4
  // Full vector field:  x_displacement * iHat + y_displacement * jHat


  // ----start_parameters_passing-----------
  char *tmp[3];
  tmp[0] = argv[0];
  tmp[1] = fe_name;
  tmp[2] = quad_name;

  if (argc == 1)  // Default parameters
    {
      argv = tmp;
    }
  else            // User-passed parameters
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
  // ----end_parameters_passing------------


  try
    {
      //=======================================
      Beam::ElasticProblem<2> elastic_problem_2d(argv[1], argv[2], degree, n_cycles_down, n_cycles_up);
      elastic_problem_2d.run();
      //=======================================
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------" << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------" << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------" << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------" << std::endl;
      return 1;
    }

  std::cout << " \n" << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= CODE ENDED CORRECTLY =====================" << std::endl;
  std::cout << "============================================================\n" << std::endl;
 
  return 0;
}

