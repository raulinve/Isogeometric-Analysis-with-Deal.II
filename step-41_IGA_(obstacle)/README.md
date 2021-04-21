# ℹ️ Step-41 IGA code (obstacle)

https://github.com/raulinve/Isogeometric-Analysis-with-Deal.II/tree/master/step-41_IGA_(obstacle)  

<br/>  

<img src="https://github.com/raulinve/Isogeometric-Analysis-with-Deal.II/blob/master/step-41_IGA_(obstacle)/doc/IMG_step-41_IGA_t6.png" alt="Obstacle example result" width="480" height="480">  

**Img. 2**: Result plot of the *obstacle* code.  

<br/>  

The following are the steps executed by the **main** function of the *obstacle* code:

<br/>  

---
### 01. Initial declarations:  
1. The default variables are:  
	```cpp
	unsigned int n_cycle       = 5;
	unsigned int degree        = 0;
	bool         h_refinement  = true;
	bool         p_refinement  = false;
	bool         k_refinement  = false;
	std::string  continuity    = "C1";
	```  

2. After the declaration of the variables, it is necessary to initialize MPI which 
    is a library for parallelization:  
	```cpp
	Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
	```  
	Note: The call to initialize MPI exists because the Trilinos library upon 
	which we build our linear solvers in this program requires it.  

3. Then, it is also needed to declare a `ConvergenceTable` object 
    and a `subdivisions` vector:  
	```cpp
	std::vector<unsigned int> subdivisions(1);
	ConvergenceTable          convergence_table;
	```

4. After these instructions, a **for loop** is started in order to perform 
**consecutive refinements** according to the selected *h-p-k* refinement type.  


<br/>  

## ↘️ for loop: start
The loop over the refining cycles starts.  
```cpp
for (unsigned int cycle=1; cycle<=n_cycle; ++cycle)
{
```

<br/>  

---
### 02. Initialization of the knot matrix and of the mults matrix:

1. Considering the refinement chosen, the following parameters are set:  

	Refinement   | `subdivisions` | `degree`
	------------ | :-----------:  | :-----------:
	h-ref        | 2^{`cycle`}    | 3
	p-ref        | 100            | `cycle`
	k-ref        | 100            | 4

2. The **knot matrix** is initialized:  
    It is declared as a vector-of-vector matrix of dimensions `1` x `subdivisions + 1`.  
    Then it is initialized as:  

    <!--- K_{0,i}={\frac{2{\cdot}i}{subdivisions}}-1  -->  

    ![equation](http://www.sciweavers.org/tex2img.php?eq=K_%7B0%2Ci%7D%3D%7B%5Cfrac%7B2%7B%5Ccdot%7Di%7D%7Bsubdivisions%7D%7D-1&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)  
    
	```cpp
	std::vector<std::vector<double> > knots(1, std::vector<double>(subdivisions[0]+1));
	for (unsigned int i=0; i<subdivisions[0]+1; ++i)
		knots[0][i] = -1+i*2./subdivisions[0];
	```  

3. The **mults matrix** is initialized:  
    It is declared as a vector-of-vector matrix of dimensions `1` x `subdivisions + 1`.  
    And it is directly initialized with all ones.  
    ```cpp
    std::vector<std::vector<unsigned int> > mults(1, std::vector<unsigned int>(subdivisions[0]+1, 1));
    ```
    Considering the selected continuity (C0-C1-C2) the matrix elements are modified accordingly:

	Continuity   | `mults[0][:]`
	------------ | -------------
	C0           | `degree`
	C1           | `degree`-1
	C2           | `degree`-2 

    If the refinement method is the **k-ref**, then the `mults` matrix is modified at each cycle as:
	```cpp
    mults[0][:] = cycle+1;
    ```

4. The ends of **open knots vectors** are enforced:  
	```cpp
	mults[0][0] = degree+1;
	mults[0][subdivisions[0]] = degree+1;
	mults.push_back(mults[0]);
	knots.push_back(knots[0]);
	```

<br/>  

---
### 03. Initialization of the IgaHandler:
After the initializations, an **IgaHandler** object is initialized:  
```cpp
IgaHandler<2,2> iga_hd2(knots, mults, degree);
```

The IgaHandler constructor template is the following:  
```cpp
template <int dim, int spacedim>
IgaHandler<dim,spacedim>::IgaHandler (const std::vector<std::vector<double> > &        knot_vectors,
                                      const std::vector<std::vector<unsigned int> > &  mults,
                                      const unsigned int                               degree)
  :
  degree            (degree),
  fe                (degree),
  dh                (tria),
  fe_sys            (FE_Bernstein<dim,spacedim>(degree), spacedim),
  dh_sys            (tria),
  map_fe            (NULL),
  local_dof_indices (fe.dofs_per_cell),
  knot_vectors      (knot_vectors),
  mults             (mults),
  local_quad        (QIterated<dim>(QTrapez<1>(), 30))
{
    // [...] see below...
}
```

Note: Of course, this is an additional step with respect to the original step-41 example.  

<br/>  

#### The constructor, when called, do the following things:  
1. Initializes all the variables as reported in the template;  

2. Initializes the `local_b_extractors`:  
	```cpp
	std::vector<std::vector<FullMatrix<double> > > local_b_extractors;  // [in the header file]

	local_b_extractors.resize(dim);
	for (int i=0; i<dim; ++i)
	    local_b_extractors[i].resize(knot_vectors[i].size()-1, IdentityMatrix(degree+1));
	```

3. Performs a **meshing of the geometry** by calling the 
`subdivided_hyper_rectangle` function:  
	```cpp
	GridGenerator::subdivided_hyper_rectangle(tria, knot_vectors);
	```  
	Whose template is:  
	```cpp
	template <int dim, int spacedim>
	void subdivided_hyper_rectangle(Triangulation<dim,spacedim>              & tria,
	                                const std::vector< std::vector<double> > & p     )
	```  
	Note: This *global function* is NOT part of any class, and its 
	implementation is contained the files: 
	"*./include/grid_generator.h*" and "*./source/grid_generator.cc*".  

	After this, the mesh is saved to a "*mesh_dim_spacedim.msh*" file.  

4. Initializes the **map_fe** object of type `MappingFEField`:  
	```cpp
	map_fe = new MappingFEField<dim,spacedim>(dh_sys, euler, mask);
	```
5. It runs the **local bezier extractor** over the knots with repetitions:  
    ```cpp
    compute_local_b_extractors(local_b_extractors[i],
                               knots_with_repetition[i]);
    ```
6. Finally it assembles the **global extractor**:  
    ```cpp
    assemble_global_extractor();
    ```
    Note: This command calls another method of the *IgaHandler* class:
    ```cpp
    template <int dim, int spacedim>
    void IgaHandler<dim,spacedim>::assemble_global_extractor()
    { /* [...] */ }
    ```

<br/>  

---
### 04. Initialization of the problem:
Then the problem is initialized by calling the **constructor** of the problem class:
```cpp
ObstacleProblem<2> obstacle_problem(iga_hd2, convergence_table);
```
Its template is:
```cpp
template <int dim>
ObstacleProblem<dim>::ObstacleProblem (IgaHandler<dim,dim> & iga_handler,
                                       ConvergenceTable &    convergence_table )

```
Note: In the original step-41 example, the method takes no parameters, indeed it was called as:
```cpp
ObstacleProblem<2> obstacle_problem;
```

<br/>  

#### The constructor only do the following things:  
1. Initializes all the variables using *iga_handler*;  
    ```cpp
    :
    iga_handler   (iga_handler),
    degree        (iga_handler.degree),
    triangulation (iga_handler.tria),
    fe            (iga_handler.fe),
    dof_handler   (iga_handler.dh),
    mappingfe     (iga_handler.map_fe),
    ```
2.  Initializes the convergence table;  
    ```cpp
    convergence_table(convergence_table)
    ```


<br/>  

---
### 05. Solution of the problem:
The problem is solved by calling the **run** method:  
```cpp
obstacle_problem.run(cycle);
```

Note: In the original step-41 example, the method takes no arguments.  
Here it is needed to pass the `cycle` parameter uniquely to compute the error using `compute_error(cycle);`.  

<br/>  

#### The run method performs the following actions:  

1. `make_grid();`  
    This method should handle the grid creation.  

    BUT: In this implementation the method is not implemented, 
    since is automatically performed inside `iga_handler`.  

    Note: In the original step-41 example the grid is produced with 
    the use of *hyper_cube*, with a fixed refinement of `refine_global(7)`.
    ```cpp
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(7);
    ```  

2. `setup_system();` (+ `assemble_mass_matrix_diagonal`)  
    This method enumerates all the degrees of freedom and sets up matrix and 
    vector objects to hold the system data.  
    
    For this purpose, *sparse matrices* (subdivided in values and pattern structures) are used.  

    The **biggest difference** between this code and the original step-41 
    is the replacement of the original `dof_handler` with the new `iga_handler`.  
    With this change, also the number of degrees of freedom changes. This information is printed on screen using:
    ```cpp
    std::cout << "   (STD dofs: " << dof_handler.n_dofs()  << ")" << std::endl
              << "   (IGA dofs: " << iga_handler.n_bspline << ")" << std::endl;
    ```  
    
    Note: This method make use of another private method in order to assemble the mass matrix:
    ```cpp
    assemble_mass_matrix_diagonal(bspline_mass_matrix);
    ```  
    This method has some differences from the original step-41, that are:  
    - In the initialization of the `FEValues<dim> fe_values` it is introduced
    a `*mappingfe,` parameter, which comes from the method `iga_handler.map_fe`.  
    - The step-41 method used to distribute the constraints from local to 
    global `constraints.distribute_local_to_global` is replaced with the corresponding `iga_handler.distribute_local_to_global`.  


<br/>  

↘️ At this point, a **for loop** starts in order to run *Newton iterations*:
```cpp
for (unsigned int iteration = 0; iteration <= solution.size(); ++iteration) {
```

<br/>  

3. `assemble_system();`  
    This method assemble the matrices and the vector producing the system to solve.  
    The implementation is very similar between the *obstacle* code and the original step-41.  
    The main differences are:  
    - In the initialization of the `FEValues<dim> fe_values` it is introduced
    a `*mappingfe,` parameter, which comes from the method `iga_handler.map_fe`.  
    - The step-41 method used to distribute the constraints from local to 
    global `constraints.distribute_local_to_global` is replaced with the corresponding `iga_handler.distribute_local_to_global`.  

    Note: See the documentation for more details.  

4. `solve();`  
    This method solve the matrix system using the preconditioned Conjugate Gradients (CG) 
    method that can be used to solve linear systems with a symmetric positive definite matrix.  
    
    Note: In the context of a *Newton method*, we are not typically interested in very high accuracy
    (because the solution of this linear problem only gives us an approximation 
    of the solution of the true nonlinear problem),
    and so we use the `ReductionControl` class that stops iterations when
    either an absolute tolerance is reached or when the residual is reduced by a certain factor.  

    
    ```cpp
    ReductionControl          reduction_control (1000, 1e-12, 1e-5);
    SolverCG<Vector<double>>  solver (reduction_control);
    precondition.initialize (bspline_system_matrix);
    solver.solve (bspline_system_matrix, bspline_solution, bspline_system_rhs, precondition);
    ```  
    Note: In the original step-41 code the solver controls were `reduction_control(100, 1e-12, 1e-3)`.  

5. `update_solution_and_constraints();`  
    This method updates the active set of constrained degrees of freedom 
    and computes an AffineConstraints object from it that can then
    be used to eliminate constrained degrees of freedom from the solution of
    the next iteration.  
    At the same time we set the constrained degrees of freedom of the 
    solution to the correct value, namely the height of the obstacle.  
    Note: In a sense, this is the central function of this example.  

    Note: In this method a lot of different small modifications were needed in 
    order to correctly treat isogeometric geometries.  
    The **most relevant modifications** are:  
    - All the `dof_handler` objects are substituted with the `iga_handler` ones;  
    - A new method is required in order to correctly match the dofs to the support points:  
        ```cpp
        iga_handler.map_dofs_to_support_points(bspline_support_points);` 
        ```  
    - A new method is required in order to project the boundary values:  
        ```cpp
        iga_handler.project_boundary_values(dirichlet_boundary,
                                            boundary_quad,
                                            bspline_constraints);
        ```  

6. `output_results(iteration);`  
    This method produces the output drawings in .vtk format.  

    Note: In order to correctly print the isogeoemetric geometries, 
    a new method is implemented able to transform iga-vectors into standard fe_spaces:  
	```cpp
	iga_handler.transform_vector_into_fe_space(
		  bspline_sol_dh,
		  bspline_solution_vector);
	```  

7. `compute_error (cycle);`  
    This method is an additon with respect to the step-41 code.  
    It is used to compute the L2-, H1- and L_Infinity-norm errors.  

<br/>  

➡️ Here the **for loop** of the *Newton iteration* ends:
```cpp
}
```

<br/>  

## ➡️ for loop: end
After the **run** method, also the loop of the main function over the refining cycles ends.  

<br/>  


---
### 06. Final procedures:
Just before the code end, the global function `print_table` is 
called to setup and save to file a *convergence table*.  
```cpp
print_table(convergence_table, h_refinement, p_refinement, k_refinement);
```



