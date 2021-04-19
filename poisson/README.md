# ℹ️ Poisson code

https://github.com/raulinve/Isogeometric-Analysis-with-Deal.II/tree/master/poisson  

<br/>  

<img src="https://github.com/raulinve/Isogeometric-Analysis-with-Deal.II/blob/master/poisson/doc/IMG_poisson_t5.png" alt="Poisson example result" width="480" height="480">  

**Img. 1**: Result plot of the *poisson* code.  

<br/>  

The following are the steps executed by the **main** function of the *poisson* code:

<br/>  

---
### 01. Definition of the main variables:  
- The variables can be passed as program arguments;  
- If no arguments are passed, the default ones are:  
  ```cpp
  char         fe_name[]     = "bernstein";
  char         quad_name[]   = "legendre";
  unsigned int degree        = 1;
  unsigned int n_cycles_down = 0;
  unsigned int n_cycles_up   = 5;
  ```

<br/>  

---
### 02. Initialization of the problem:
The problem is initialized by calling the **constructor** of the problem class:
```cpp
Laplace<2> laplace_problem_2d(argv[1], argv[2], degree, n_cycles_down, n_cycles_up);
```
Its template is:
```cpp
template <int dim>
Laplace<dim>::Laplace (const std::string   fe_name,
                       const std::string   quadrature_name,
                       const unsigned int  degree,
                       const unsigned int  n_cycles_low,
                       const unsigned int  n_cycles_up );
```

<br/>  

#### The constructor do the following things:  
1. Initializes all the passed variables;  
2. Specifies the polynomial degree of the finite elements (in this case *NULL*);  
    `fe(NULL),`  
3. Associates the DoFHandler to the triangulation;  
    `dof_handler (triangulation)`  

Note: Until here the code is identical to the original step-4 example.  

4. Then the constructor initializes the **quadrature formula**:  
    ```cpp
    if (quadrature_name == "legendre")  {    // <-[default]
      matrix_quad   = QGauss<dim>  (degree+1);
      boundary_quad = QGauss<dim-1>(degree+2);
      error_quad    = QGauss<dim>  (degree+3);
    }
    else if (quadrature_name == "lobatto")  {
      matrix_quad   = QGaussLobatto<dim>  (degree+2);
      boundary_quad = QGaussLobatto<dim-1>(degree+3);
      error_quad    = QGaussLobatto<dim>  (degree+4);
    }
    ```

5. And it also initializes the **FiniteElement** type:  
    ```cpp
    if (fe_name == "bernstein")
      fe = new FE_Bernstein<dim>(degree);
    else if (fe_name == "lagrange")    // <-[default]
      fe = new FE_Q<dim>(degree);
    else if (fe_name == "lobatto")
      fe = new FE_Q<dim>(QGaussLobatto<1>(degree+1));
    ```
    Note: This initialization at runtime requires to declare the finite element object as a pointer, since it is not known yet at compile time and it does not support assignement operators.  


<br/>  

---
### 03. Solution of the problem:
The problem is solved by calling the **run** method:  
```cpp
laplace_problem_2d.run ();
```

<br/>  

#### The run method performs the following actions:  

1. `make_grid();`  
    This method handle the grid creation.  
    Note: This method is invoked only once, when the first cycle runs.  
    The grid is produced with the use of *hyper_cube*:  
    ```cpp
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(n_cycles_low+1);
    ```
    In the original step-4 example the refinement is fixed at `refine_global(4)`.  

2. `refine_grid();`  
    This method is an additon with respect to the step-4 code.  
    It is used to refine the grid at each cycle using the library command:  
    ```cpp
    triangulation.refine_global (1);
    ```

3. `setup_system();`  
    This method enumerates all the degrees of freedom and sets up matrix and 
    vector objects to hold the system data.  
    For this purpose, a sparse matrix (subdivided in values and pattern structures) is used.  
    The only difference between this code and the original step-4 is the use of `*fe` instead of `fe`.  

4. `assemble_system();`  
    This method assemble the matrices and the vector producing the system to solve.  
    The implementation is very similar between the *poisson* code and the original step-41.  
    Note: See the documentation for more details.  

5. `solve();`  
    This method solve the matrix system using the preconditioned Conjugate Gradients (CG) 
    method that can be used to solve linear systems with a symmetric positive definite matrix.  
    ```cpp
    SolverControl            solver_control(100000, 1e-14);
    SolverCG<Vector<double>> solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
    ```  
    Note: In the original step-4 code the solver controls were `(1000, 1e-12)`.  

6. `output_results(cycle);`  
    This method produces the output drawings in .vtk format.  

7. `process_solution(cycle);`  
    This method is an additon with respect to the step-4 code.  
    It is used to compute the L2- and H1-norm errors.  


8. `print_table(cycle);`  
    This method is an additon with respect to the step-4 code.  
    It is used to setup and save to file a convergence table.  





