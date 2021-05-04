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

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>   // [ ADD ]
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// In this example, we need vector-valued finite elements. The support for
// these can be found in the following include file:
#include <deal.II/fe/fe_system.h>
// We will compose the vector-valued finite elements from regular Q1 elements
// which can be found here, as usual:
#include <deal.II/fe/fe_q.h>

#include <deal.II/fe/fe.h>              // IGA
#include <deal.II/fe/fe_nothing.h>      // IGA
#include <deal.II/fe/fe_bernstein.h>    // IGA
#include <deal.II/lac/sparse_ilu.h>     // IGA

// This again is C++:
#include <fstream>
#include <iostream>

#include <filesystem>    // IGA

// The last step is as in previous programs. In particular, just like in
// step-7, we pack everything that's specific to this program into a namespace
// of its own.
using namespace dealii;
namespace Step8
{
  //using namespace dealii;

  bool original_problem_8 = true;	// Default: "true" [SHOULD BE REMOVED NEXT]

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
                        const unsigned int  component = 0) const 
  {
    return 0*p(0)*component;
  }
};


//====================================================
  // @sect3{Right hand side values}

  // Before going over to the implementation of the main class, we declare and
  // define the function which describes the right hand side. This time, the
  // right hand side is vector-valued, as is the solution, so we will describe
  // the changes required for this in some more detail.
  //
  // To prevent cases where the return vector has not previously been set to
  // the right size we test for this case and otherwise throw an exception at
  // the beginning of the function. Note that enforcing that output arguments
  // already have the correct size is a convention in deal.II, and enforced
  // almost everywhere. The reason is that we would otherwise have to check at
  // the beginning of the function and possibly change the size of the output
  // vector. This is expensive, and would almost always be unnecessary (the
  // first call to the function would set the vector to the right size, and
  // subsequent calls would only have to do redundant checks). In addition,
  // checking and possibly resizing the vector is an operation that can not be
  // removed if we can't rely on the assumption that the vector already has
  // the correct size; this is in contract to the Assert call that is
  // completely removed if the program is compiled in optimized mode.
  //
  // Likewise, if by some accident someone tried to compile and run the
  // program in only one space dimension (in which the elastic equations do
  // not make much sense since they reduce to the ordinary Laplace equation),
  // we terminate the program in the second assertion. The program will work
  // just fine in 3d, however.
  template <int dim>
  void right_hand_side(const std::vector<Point<dim>> &points,
                       std::vector<Tensor<1, dim>> &  values)
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    Assert(dim >= 2, ExcNotImplemented());

	if(original_problem_8) {
		// The rest of the function implements computing force values. We will use
		// a constant (unit) force in x-direction located in two little circles
		// (or spheres, in 3d) around points (0.5,0) and (-0.5,0), and y-force in
		// an area around the origin; in 3d, the z-component of these centers is
		// zero as well.
		//
		// For this, let us first define two objects that denote the centers of
		// these areas. Note that upon construction of the Point objects, all
		// components are set to zero.
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
		// Points of application:
		Point<dim> point_1, point_2;
		point_1(0) =  46;   point_1(1) =  50;
		//point_2(0) = -0.5;   point_2(1) = -0.5;
	
		double radius = 2;    // Radius of the force footprint
		double force  = -0.1;     // Absolute modulus of the forces

		for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
		  { 
		    if (((points[point_n] - point_1).norm_square() < radius * radius))
		      values[point_n][1] = force;
		    else
		      values[point_n][1] = 0.0;
            /*
		    // If <code>points[point_n]</code> is in a circle (sphere) of radius
		    // 0.2 around one of these points, then set the force in x-direction
		    // to one, otherwise to zero:
		    if (((points[point_n] - point_1).norm_square() < radius * radius) ||
		        ((points[point_n] - point_2).norm_square() < radius * radius))
		      values[point_n][0] = force;
		    else
		      values[point_n][0] = 0.0;

		    // Likewise, if <code>points[point_n]</code> is in the vicinity of the
		    // origin, then set the y-force to one, otherwise to zero:
		    if (points[point_n].norm_square() < radius * radius)
		      values[point_n][1] = force;
		    else
		      values[point_n][1] = 0.0;
            */
		  }
	}
  }


//====================================================
  // @sect3{The <code>ElasticProblem</code> class template}

  // The main class is, except for its name, almost unchanged with respect to
  // the step-6 example.
  //
  // The only change is the use of a different class for the <code>fe</code>
  // variable: Instead of a concrete finite element class such as FE_Q, we now
  // use a more generic one, FESystem. In fact, FESystem is not really a
  // finite element itself in that it does not implement shape functions of
  // its own. Rather, it is a class that can be used to stack several other
  // elements together to form one vector-valued finite element. In our case,
  // we will compose the vector-valued element of <code>FE_Q(1)</code>
  // objects, as shown below in the constructor of this class.
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
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void refine_grid();
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

    Quadrature<dim>      matrix_quad;        // IGA
    Quadrature<dim>      error_quad;         // IGA
    Quadrature<dim-1>    boundary_quad;      // IGA
  };



  // @sect3{The <code>ElasticProblem</code> class implementation}

  // @sect4{ElasticProblem::ElasticProblem constructor}

  // Following is the constructor of the main class. As said before, we would
  // like to construct a vector-valued finite element that is composed of
  // several scalar finite elements (i.e., we want to build the vector-valued
  // element so that each of its vector components consists of the shape
  // functions of a scalar element). Of course, the number of scalar finite
  // elements we would like to stack together equals the number of components
  // the solution function has, which is <code>dim</code> since we consider
  // displacement in each space direction. The FESystem class can handle this:
  // we pass it the finite element of which we would like to compose the
  // system of, and how often it shall be repeated:

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
    fe              (NULL),        //fe(FE_Q<dim>(1), dim)
    dof_handler     (triangulation)
  {
    // IGA: This section is entirely new.
    if (quadrature_name == "legendre")     // [DEFAULT]
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
    else if (fe_name == "lagrange")        // [DEFAULT]
      fe = new FESystem<dim>(FE_Q<dim>(degree), dim);
    else if (fe_name == "lobatto")
      fe = new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(degree+1)), dim);
    else
	  AssertThrow(false, ExcMessage("FE not supported"));

  }
  // In fact, the FESystem class has several more constructors which can
  // perform more complex operations than just stacking together several
  // scalar finite elements of the same type into one; we will get to know
  // these possibilities in later examples.

/// Main problem class: Destructor
	template <int dim>
	ElasticProblem<dim>::~ElasticProblem ()
	{
	  dof_handler.clear();
	  if (fe)
		delete fe;
	}


// On to the first of the private member functions. Here we create the
// triangulation of the domain, for which we choose a scaled an anisotripically
// discretised rectangle which is subsequently transformed into the correct
// of the Cook cantilever. Each relevant boundary face is then given a boundary
// ID number.
//
// We then determine the volume of the reference configuration and print it
// for comparison.

  template <int dim>
  Point<dim> grid_y_transform (const Point<dim> &pt_in)                             // [ ADD ]
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

  template <int dim>
  void ElasticProblem<dim>::make_grid()
  {
	if(original_problem_8) {
	    GridGenerator::hyper_cube(triangulation, -1, 1);
	    //triangulation.refine_global(4);
        triangulation.refine_global(n_cycles_low+4);    // IGA
	}
	else {
		// Divide the beam, but only along the x- and y-coordinate directions
		unsigned int elements_per_edge = 16;
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
		          cell->face(face)->set_boundary_id(1); // -X faces
		        else if (std::abs(cell->face(face)->center()[0] - 48.0) < tol_boundary)
		          cell->face(face)->set_boundary_id(11); // +X faces
		        else if (std::abs(std::abs(cell->face(face)->center()[0]) - 0.5) < tol_boundary)
		          cell->face(face)->set_boundary_id(2); // +Z and -Z faces
		      }

		// Transform the hyper-rectangle into the beam shape
		GridTools::transform(&grid_y_transform<dim>, triangulation);


		//GridTools::scale(1e-3, triangulation);   // parameters.scale

		/*
		vol_reference = GridTools::volume(triangulation);
		vol_current = vol_reference;
		std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;
		*/
	} 
  }


  // @sect4{ElasticProblem::setup_system}

  // Setting up the system of equations is identical to the function used in
  // the step-6 example. The DoFHandler class and all other classes used here
  // are fully aware that the finite element we want to use is vector-valued,
  // and take care of the vector-valuedness of the finite element
  // themselves. (In fact, they do not, but this does not need to bother you:
  // since they only need to know how many degrees of freedom there are per
  // vertex, line and cell, and they do not ask what they represent,
  // i.e. whether the finite element under consideration is vector-valued or
  // whether it is, for example, a scalar Hermite element with several degrees
  // of freedom on each vertex).
  template <int dim>
  void ElasticProblem<dim>::setup_system()
  {

    dof_handler.distribute_dofs(*fe);        // IGA: "*"
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

	// IMPOSE CONSTRAINTS: -------------
    constraints.clear();

	if(original_problem_8) {
      std::cout << "A" << std::endl; //XXX [DEBUG FLAG]

      std::map< types::boundary_id, const Function< dim > *>  dirichlet_boundary;
      BoundaryValues<dim> boundary_funct;                                 // IGA
      dirichlet_boundary[0] = &boundary_funct;                            // IGA

	  //std::map<types::global_dof_index,double> boundary_values;           // IGA

      DoFTools::make_hanging_node_constraints(dof_handler, constraints);

      std::cout << "B" << std::endl; //XXX [DEBUG FLAG]

      // NOTE: "interpolate" method works fine with std Lagrange polynomial, BUT FAILS with Bernstein. 
      /*VectorTools::interpolate_boundary_values (dof_handler,       
                                                0,
                                                Functions::ZeroFunction<dim>(dim),
                                                constraints);*/

      // NOTE: This is why it is needed to replece it with the following "project" method.
      VectorTools::project_boundary_values (dof_handler,         // IGA: from "VectorTools::interpolate_boundary_values(..."
                                            dirichlet_boundary,  // IGA: from "0"
                                            boundary_quad,       // IGA: from "Functions::ZeroFunction<dim>(dim)"
                                            constraints);        // IGA: from "constraints"  // TRY: "boundary_values"

      std::cout << "C" << std::endl; //XXX [DEBUG FLAG]
	}
	else {
	  // Fixed left hand side of the beam
      const int boundary_id = 1;
	  const FEValuesExtractors::Vector u_fe(0);
      VectorTools::interpolate_boundary_values(dof_handler,  // dof_handler_ref
                                               boundary_id,
                                               ZeroFunction<dim>(dim),
                                               constraints,
                                               fe->component_mask(u_fe));    // IGA: from "fe.component_mask(u_fe));"
	}
    constraints.close();
	// ---------------------------------

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);

    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    //solution.reinit (dof_handler.n_dofs());    // IGA
    //system_rhs.reinit (dof_handler.n_dofs());  // IGA

  }


  // @sect4{ElasticProblem::assemble_system}

  // The big changes in this program are in the creation of matrix and right
  // hand side, since they are problem-dependent. We will go through that
  // process \step-by-step, since it is a bit more complicated than in previous
  // examples.
  //
  // The first parts of this function are the same as before, however: setting
  // up a suitable quadrature formula, initializing an FEValues object for the
  // (vector-valued) finite element we use as well as the quadrature object,
  // and declaring a number of auxiliary arrays. In addition, we declare the
  // ever same two abbreviations: <code>n_q_points</code> and
  // <code>dofs_per_cell</code>. The number of degrees of freedom per cell we
  // now obviously ask from the composed finite element rather than from the
  // underlying scalar Q1 element. Here, it is <code>dim</code> times the
  // number of degrees of freedom per cell of the Q1 element, though this is
  // not explicit knowledge we need to care about:
  template <int dim>
  void ElasticProblem<dim>::assemble_system()
  {
    std::cout << "   > ASSEMBLING THE SYSTEM (wait) ... " << std::endl;
    // QGauss<dim> quadrature_formula(fe.degree + 1);    // IGA deleted

    FEValues<dim> fe_values(*fe, matrix_quad,     // IGA: "*", "quadrature_formula" to "matrix_quad"
                            update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe->dofs_per_cell;   // IGA: "->"
    const unsigned int n_q_points    = matrix_quad.size();  // IGA from "quadrature_formula.size();"

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Addition in case variable coefficients are needed:
    std::vector<double> lambda_values(n_q_points);        // [not needed]
    std::vector<double> mu_values(n_q_points);            // [not needed]
    Functions::ConstantFunction<dim> lambda(1.), mu(1.);  // [not needed]

    // Like the two constant functions above, we will call the function
    // right_hand_side just once per cell to make things simpler.
    std::vector<Tensor<1, dim>> rhs_values(n_q_points);


    typename DoFHandler<dim>::active_cell_iterator    // IGA new
    cell = dof_handler.begin_active(),                // IGA new
    endc = dof_handler.end();                         // IGA new
    // Now we can begin with the loop over all cells:
    for (; cell!=endc; ++cell)    // IGA (from:  "for (const auto &cell : dof_handler.active_cell_iterators())")
      {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs    = 0;

        // Next we get the values of the coefficients at the quadrature points. 
        lambda.value_list(fe_values.get_quadrature_points(), lambda_values);    // [not needed]
        mu.value_list(fe_values.get_quadrature_points(), mu_values);            // [not needed]
        right_hand_side(fe_values.get_quadrature_points(), rhs_values);         // [not needed]

        // Then assemble the entries of the local stiffness matrix and right
        // hand side vector. This follows almost one-to-one the pattern
        // described in the introduction of this example.  One of the few
        // comments in place is that we can compute the number
        // <code>comp(i)</code>, i.e. the index of the only nonzero vector
        // component of shape function <code>i</code> using the
        // <code>fe.system_to_component_index(i).first</code> function call
        // below.
        //
        // (By accessing the <code>first</code> variable of the return value
        // of the <code>system_to_component_index</code> function, you might
        // already have guessed that there is more in it. In fact, the
        // function returns a <code>std::pair@<unsigned int, unsigned
        // int@></code>, of which the first element is <code>comp(i)</code>
        // and the second is the value <code>base(i)</code> also noted in the
        // introduction, i.e.  the index of this shape function within all the
        // shape functions that are nonzero in this component,
        // i.e. <code>base(i)</code> in the diction of the introduction. This
        // is not a number that we are usually interested in, however.)              // [FROM HERE !!!]
        //
        // With this knowledge, we can assemble the local matrix
        // contributions:
        for (unsigned int i=0; i<dofs_per_cell; ++i)  // IGA from: "const unsigned int i : fe_values.dof_indices()"  
          {
            const unsigned int component_i = fe->system_to_component_index(i).first;    // IGA: "->"

            for (unsigned int j=0; j<dofs_per_cell; ++j)  // IGA from: "const unsigned int j : fe_values.dof_indices()"
              {
                const unsigned int component_j =fe->system_to_component_index(j).first;    // IGA: "->"

                for (unsigned int q_point=0; q_point<n_q_points; ++q_point)    // IGA from: "const unsigned int q_point : fe_values.quadrature_point_indices())"
                  {
                    cell_matrix(i, j) +=
                      // The first term is $\lambda \partial_i u_i, \partial_j
                      // v_j) + (\mu \partial_i u_j, \partial_j v_i)$. Note
                      // that <code>shape_grad(i,q_point)</code> returns the
                      // gradient of the only nonzero component of the i-th
                      // shape function at quadrature point q_point. The
                      // component <code>comp(i)</code> of the gradient, which
                      // is the derivative of this only nonzero vector
                      // component of the i-th shape function with respect to
                      // the comp(i)th coordinate is accessed by the appended
                      // brackets.
                      (                                                  //
                        (fe_values.shape_grad(i, q_point)[component_i] * //
                         fe_values.shape_grad(j, q_point)[component_j] * //
                         lambda_values[q_point])                         //
                        +                                                //
                        (fe_values.shape_grad(i, q_point)[component_j] * //
                         fe_values.shape_grad(j, q_point)[component_i] * //
                         mu_values[q_point])                             //
                        +                                                //
                        // The second term is $(\mu \nabla u_i, \nabla
                        // v_j)$. We need not access a specific component of
                        // the gradient, since we only have to compute the
                        // scalar product of the two gradients, of which an
                        // overloaded version of <tt>operator*</tt> takes
                        // care, as in previous examples.
                        //
                        // Note that by using the <tt>?:</tt> operator, we only
                        // do this if <tt>component_i</tt> equals
                        // <tt>component_j</tt>, otherwise a zero is added
                        // (which will be optimized away by the compiler).
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
        for (unsigned int i=0; i<dofs_per_cell; ++i)    // IGA: from "const unsigned int i : fe_values.dof_indices())"
          {
            const unsigned int component_i = fe->system_to_component_index(i).first;    // IGA: "->"

            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)    // IGA from: "const unsigned int q_point : fe_values.quadrature_point_indices())"
              cell_rhs(i) += fe_values.shape_value(i, q_point) *
                             rhs_values[q_point][component_i] *
                             fe_values.JxW(q_point);
          }

        // The transfer from local degrees of freedom into the global matrix
        // and right hand side vector does not depend on the equation under
        // consideration, and is thus the same as in all previous
        // examples.
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }



  // @sect4{ElasticProblem::solve}

  // The solver does not care about where the system of equations comes, as
  // long as it stays positive definite and symmetric (which are the
  // requirements for the use of the CG solver), which the system indeed
  // is. Therefore, we need not change anything.
  template <int dim>
  void ElasticProblem<dim>::solve()
  {
    SolverControl            solver_control(1000, 1e-12);
    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    std::cout << "   Memory consumption " << system_matrix.memory_consumption()
              << " bytes" << std::endl;

    std::cout << "   > SOLVING THE SYSTEM (wait) ... " << std::endl;
    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    std::cout << "   " << solver_control.last_step()
              << " CG iterations needed to obtain convergence." << std::endl;

    cg_iter = solver_control.last_step();
  }


  // @sect4{ElasticProblem::refine_grid}

  // The function that does the refinement of the grid is the same as in the
  // step-6 example. The quadrature formula is adapted to the linear elements
  // again. Note that the error estimator by default adds up the estimated
  // obtained from all components of the finite element solution, i.e., it
  // uses the displacement in all directions with the same weight. If we would
  // like the grid to be adapted to the x-displacement only, we could pass the
  // function an additional parameter which tells it to do so and do not
  // consider the displacements in all other directions for the error
  // indicators. However, for the current problem, it seems appropriate to
  // consider all displacement components with equal weight.
  template <int dim>
  void ElasticProblem<dim>::refine_grid()
  {
    std::cout << "   Refining grid." << std::endl;
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate(dof_handler,
                                       QGauss<dim - 1>(fe.degree + 1),
                                       {},
                                       solution,
                                       estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);

    triangulation.execute_coarsening_and_refinement();
  }


  // @sect4{ElasticProblem::output_results}

  // The output happens mostly as has been shown in previous examples
  // already. The only difference is that the solution function is vector
  // valued. The DataOut class takes care of this automatically, but we have
  // to give each component of the solution vector a different name.
  //
  // To do this, the DataOut::add_vector() function wants a vector of
  // strings. Since the number of components is the same as the number
  // of dimensions we are working in, we use the <code>switch</code>
  // statement below.
  //
  // We note that some graphics programs have restriction on what
  // characters are allowed in the names of variables. deal.II therefore
  // supports only the minimal subset of these characters that is supported
  // by all programs. Basically, these are letters, numbers, underscores,
  // and some other characters, but in particular no whitespace and
  // minus/hyphen. The library will throw an exception otherwise, at least
  // if in debug mode.
  //
  // After listing the 1d, 2d, and 3d case, it is good style to let the
  // program die if we run upon a case which we did not consider. Remember
  // that the Assert macro generates an exception if the condition in the
  // first parameter is not satisfied. Of course, the condition
  // <code>false</code> can never be satisfied, so the program will always
  // abort whenever it gets to the default statement:
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
    //std::cout << "   Solution: " << solution.size() << " [of 578]" << std::endl;                        //XXX
    //unsigned int k = 0;                                                                                 //XXX
    //for (auto i : solution) {                                                                           //XXX
    //  std::cout << "   > " << std::setfill('0') << std::setw(3) << ++k << " : " << i << std::endl; }    //XXX
    data_out.add_data_vector(solution, solution_names);
    
    data_out.build_patches();

	std::string filename = "solution-" + std::to_string(dim) + "d-";
	filename += ('0' + cycle);
	filename += ".vtk";

    std::string relpath = "RESULTS/" + filename;    // ADD
    std::ofstream output (relpath);
    data_out.write_vtk(output);

  }



  // @sect4{ElasticProblem::run}

  // The <code>run</code> function does the same things as in step-6, for
  // example. This time, we use the square [-1,1]^d as domain, and we refine
  // it globally four times before starting the first iteration.
  //
  // The reason for refining is a bit accidental: we use the QGauss
  // quadrature formula with two points in each direction for integration of the
  // right hand side; that means that there are four quadrature points on each
  // cell (in 2D). If we only refine the initial grid once globally, then there
  // will be only four quadrature points in each direction on the
  // domain. However, the right hand side function was chosen to be rather
  // localized and in that case, by pure chance, it happens that all quadrature
  // points lie at points where the right hand side function is zero (in
  // mathematical terms, the quadrature points happen to be at points outside
  // the <i>support</i> of the right hand side function). The right hand side
  // vector computed with quadrature will then contain only zeroes (even though
  // it would of course be nonzero if we had computed the right hand side vector
  // exactly using the integral) and the solution of the system of
  // equations is the zero vector, i.e., a finite element function that is zero
  // everywhere. In a sense, we
  // should not be surprised that this is happening since we have chosen
  // an initial grid that is totally unsuitable for the problem at hand.
  //
  // The unfortunate thing is that if the discrete solution is constant, then
  // the error indicators computed by the KellyErrorEstimator class are zero
  // for each cell as well, and the call to
  // Triangulation::refine_and_coarsen_fixed_number() will not flag any cells
  // for refinement (why should it if the indicated error is zero for each
  // cell?). The grid in the next iteration will therefore consist of four
  // cells only as well, and the same problem occurs again.
  //
  // The conclusion needs to be: while of course we will not choose the
  // initial grid to be well-suited for the accurate solution of the problem,
  // we must at least choose it such that it has the chance to capture the
  // important features of the solution. In this case, it needs to be able to
  // see the right hand side. Thus, we refine globally four times. (Any larger
  // number of global refinement steps would of course also work.)
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

    // Perform the cycles:
  for (unsigned int cycle = n_cycles_low; cycle < n_cycles_up+1; ++cycle)  // IGA cycles
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;
        if (cycle == n_cycles_low)
          {
            make_grid();                      // IGA ADD (+custom geometry)
          }
        else 
          {
            if(original_problem_8) {
              //refine_grid();                    // smart refinement
              triangulation.refine_global(1);   // IGA: constant refinement
            }
            else {
              triangulation.refine_global(1);   // IGA: constant refinement
            }
          }

      std::cout << "   Number of active cells: "  //
                << triangulation.n_active_cells() //
                << std::endl                      //
                << "   Total number of cells: "   //
                << triangulation.n_cells()        //
                << std::endl;

        setup_system();
        std::cout << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl;
        assemble_system();
        solve();
        output_results(cycle);
      }
  }
} // namespace Step8

// @sect3{The <code>main</code> function}

// After closing the <code>Step8</code> namespace in the last line above, the
// following is the main function of the program and is again exactly like in
// step-6 (apart from the changed class names, of course).
int main(int argc, char **argv)
{
  std::cout << " " << std::endl;
  std::cout << "============================================================" << std::endl;
  std::cout << "================= THE CODE IS RUNNING ======================" << std::endl;
  std::cout << "============================================================" << std::endl;

  // ----start_parameters_passing-----------
  // DEFAULT VALUES:
  //--------------------------------------------------------
  char fe_name[]             = "lagrange"; // “lagrange”
  char quad_name[]           = "legendre";  // “legendre”
  unsigned int degree        = 1;           // 1
  unsigned int n_cycles_down = 0;           // 0
  unsigned int n_cycles_up   = 0;           // 2
  //--------------------------------------------------------
  // ./step-8 lagrange legendre 1 0 0
  // ./step-8 bernstein legendre 1 0 0
  // Full vector field:  x_displacement * iHat + y_displacement * jHat

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
      Step8::ElasticProblem<2> elastic_problem_2d(argv[1], argv[2], degree, n_cycles_down, n_cycles_up);
      elastic_problem_2d.run();
      //=======================================
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

