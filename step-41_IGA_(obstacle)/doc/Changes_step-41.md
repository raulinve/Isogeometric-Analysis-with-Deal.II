# Main differences between original step-41 and IGA-Obstacle

<br/>  


## MAIN:

```cpp
std::vector<unsigned int> subdivisions(1);           // IGA
for (unsigned int cycle=1; cycle<=n_cycle; ++cycle)  // IGA
  ... [IGA: build knots/mults/deg based on refinement and continuity ] ...

  IgaHandler<2,2> iga_hd2(knots, mults, degree);     // IGA
  ObstacleProblem<2> obstacle_problem(iga_hd2);      // IGA(passing iga_h2d)
  obstacle_problem.run();
```

<br/>  


## CLASS: ObstacleProblem

---
### ↘️ CONSTRUCTOR:

FROM:
```cpp
  ObstacleProblem<dim>::ObstacleProblem()
    : fe(1)
    , dof_handler(triangulation)
  {}
```
TO:
```cpp
  ObstacleProblem<dim>::ObstacleProblem (IgaHandler<dim,dim> &iga_handler)
    :
    iga_handler(iga_handler),          // IGA new
    degree(iga_handler.degree),        // IGA new
    triangulation(iga_handler.tria),   // IGA new
    fe (iga_handler.fe),               // IGA (mod)
    dof_handler (iga_handler.dh),      // IGA (mod)
    mappingfe(iga_handler.map_fe)      // IGA new
  {}
```

---
### ↘️ ADDED DESTRUCTOR:

```cpp
  template <int dim>
  ObstacleProblem<dim>::~ObstacleProblem ()
  {
	bspline_system_matrix.clear();
	bspline_complete_system_matrix.clear();
	bspline_mass_matrix.clear();
	if (mappingfe)
	  delete mappingfe;
  }
```

---
### ↘️ run ()
```cpp
// IDENTICO
```

---
### ↘️ output_results(const unsigned int iteration) const

FROM:
```cpp
TrilinosWrappers::MPI::Vector active_set_vector(dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
for (const auto index : active_set)
  active_set_vector[index] = 1.;
```
TO:
```cpp
Vector<double>  active_set_vector;  // IGA  [ MUST BE TAKEN OUT from output_results ]
```

```cpp
DataOut<dim> data_out;
data_out.attach_dof_handler(dof_handler);
  +[IGA Transformations]                                    // IGA Only
data_out.add_data_vector(solution, "displacement");
  +[IGA Transformations]                                    // IGA Only
data_out.add_data_vector (obstacle_vector, "obstacle");     // IGA Only
  +[IGA Transformations]                                    // IGA Only
data_out.add_data_vector(active_set_vector, "active_set");
  +[IGA Transformations]                                    // IGA Only
data_out.add_data_vector(contact_force, "lambda");
data_out.build_patches();
```

---
### ↘️ solve ()

FROM:
```cpp
ReductionControl        reduction_control(100, 1e-12, 1e-3);
```
TO:
```cpp
ReductionControl        reduction_control (1000, 1e-12, 1e-5);
```

```cpp
TrilinosWrappers::PreconditionAMG   precondition;      // [ MUST BE DECLARED IN MAIN CLASS !! ]
//cg_iter = reduction_control.last_step();             // [ ?? ?? ?? ]
```

---
### ↘️ update_solution_and_constraints ()

FROM:
```cpp
TrilinosWrappers::MPI::Vector lambda(complete_index_set(dof_handler.n_dofs()));
complete_system_matrix.residual(lambda, solution, complete_system_rhs);
```
TO:
```cpp
Vector<double>       bspline_lambda;     // IGA [ MUST BE TAKEN OUT from update_solution_and_constraints ]
bspline_lambda = 0;
bspline_complete_system_matrix.residual(bspline_lambda, bspline_solution, bspline_complete_system_rhs);
```

FROM:
```cpp
contact_force = lambda;
contact_force.scale(diagonal_of_mass_matrix);
contact_force *= -1;
```
TO:
```cpp
for (unsigned int i=0; i<bspline_contact_force.size(); ++i)
  bspline_contact_force[i] = bspline_lambda[i]/diagonal_of_mass_matrix[i];
bspline_contact_force *= -1;
```

```cpp
constraints.clear();
active_set.clear();
active_set_vector.reinit(0);                            // IGA Only
active_set_vector.reinit(iga_handler.n_bspline);        // IGA Only

const Obstacle<dim> obstacle;
```

FROM:
```cpp
std::vector<bool>   dof_touched(dof_handler.n_dofs(), false);
for (const auto &cell : dof_handler.active_cell_iterators())
  for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
    {
      Assert(dof_handler.get_fe().dofs_per_cell ==
               GeometryInfo<dim>::vertices_per_cell,
             ExcNotImplemented());

      const unsigned int dof_index = cell->vertex_dof_index(v, 0);
```
TO:
```cpp
std::vector<bool>   dof_touched (iga_handler.n_bspline, false);

std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);
std::vector<Point<dim> > bspline_support_points(iga_handler.n_bspline);
iga_handler.map_dofs_to_support_points(bspline_support_points);

typename DoFHandler<dim>::active_cell_iterator
cell = dof_handler.begin_active(),
endc = dof_handler.end();
for (; cell!=endc; ++cell) {
    local_dof_indices = iga_handler.iga_objects[cell].dof_indices;
    for (unsigned int v=0; v<fe.dofs_per_cell; ++v)
      {
        const unsigned int dof_index = local_dof_indices[v];
```

\[SOME SMALL DIFFERENCES]...

FROM:
```cpp
constraints.close();
```
TO:
```cpp
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
```

---
### ↘️ assemble_mass_matrix_diagonal (SparseMatrix<double> &mass_matrix)

FROM:
```cpp
Assert(fe.degree == 1, ExcNotImplemented());
const QTrapez<dim> quadrature_formula;
FEValues<dim>      fe_values(fe,
                             quadrature_formula,
                             update_values | update_JxW_values);
```
TO:
```cpp
QIterated<dim>     quadrature_formula(QTrapez<1>(),fe.degree);
FEValues<dim>      fe_values (*mappingfe,    //IGA (Add Parameter)
                              fe,
                              quadrature_formula,
                              update_values | update_JxW_values);
```
\[...]

FROM:
```cpp
 for (const auto &cell : dof_handler.active_cell_iterators())
```
TO:
```cpp
 typename DoFHandler<dim>::active_cell_iterator
 cell = dof_handler.begin_active(),
 endc = dof_handler.end();
 for (; cell!=endc; ++cell)
```
SAME:
```cpp
{
fe_values.reinit(cell);
cell_matrix = 0;
for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    for (unsigned int j=0; j<dofs_per_cell; ++j)    // IGA Only [!!]
      cell_matrix(i, i) += [SAME...]
```
FROM:
```cpp
cell->get_dof_indices(local_dof_indices);
constraints.distribute_local_to_global(cell_matrix,
                                       local_dof_indices,
                                       mass_matrix);
```
TO:
```cpp
iga_handler.distribute_local_to_global(cell_matrix,
                                   cell,
                                   mass_matrix,
                                   bspline_constraints);
}
```

---
### ↘️ assemble_system ()

FROM:
```cpp
FEValues<dim> fe_values (fe,
                         quadrature_formula,
                         update_values | update_gradients |
                         update_quadrature_points | update_JxW_values);
```
TO:
```cpp
FEValues<dim> fe_values (*mappingfe, 
                         fe,
                         quadrature_formula,
                         update_values   | update_gradients |
                         update_quadrature_points | update_JxW_values);
```

FROM:
```cpp
 for (const auto &cell : dof_handler.active_cell_iterators())
  {
    [SAME ...]
    cell->get_dof_indices(local_dof_indices);
    constraints.distribute_local_to_global(cell_matrix,
                                           cell_rhs,
                                           local_dof_indices,
                                           system_matrix,
                                           system_rhs,
                                           true);
```
TO:
```cpp
 typename DoFHandler<dim>::active_cell_iterator
 cell = dof_handler.begin_active(),
 endc = dof_handler.end();
 for (; cell!=endc; ++cell)
  {
    [SAME ...]
    iga_handler.distribute_local_to_global(cell_matrix,
                                           cell_rhs,
                                           cell,
                                           bspline_system_matrix,
                                           bspline_system_rhs,
                                           bspline_constraints);
```

---
### ↘️ setup_system()

FROM:
```cpp
active_set.set_size (dof_handler.n_dofs());
```

TO:
```cpp
active_set.set_size (iga_handler.n_bspline);
active_set_vector.reinit(iga_handler.n_bspline);
```

NOTE:
```cpp
std::cout << "   (STD dofs: " << dof_handler.n_dofs()  << ")" << std::endl
          << "   (IGA dofs: " << iga_handler.n_bspline << ")" << std::endl;
```


FROM:
```cpp
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
```

TO:
```cpp
DynamicSparsityPattern bspline_sp(iga_handler.n_bspline);    // little mod
iga_handler.make_sparsity_pattern (bspline_sp);              // little mod

SparsityPattern      sparsity_bspline;                      // IGA [declared from setup_system]
sparsity_bspline.copy_from(bspline_sp);                     // iga
bspline_system_matrix.reinit(sparsity_bspline);             // iga
bspline_complete_system_matrix.reinit(sparsity_bspline);    // iga

bspline_solution.reinit (iga_handler.n_bspline);
bspline_system_rhs.reinit (iga_handler.n_bspline);
bspline_complete_system_rhs.reinit (iga_handler.n_bspline);
bspline_contact_force.reinit (iga_handler.n_bspline);
bspline_lambda.reinit (iga_handler.n_bspline);              //IGA Only

SparseMatrix<double> bspline_mass_matrix;             // IGA [ MUST BE TAKEN OUT from setup_system ]
Vector<double>       mass_lumping;                    // IGA [ MUST BE TAKEN OUT from setup_system ]
Vector<double>       bspline_mass_lumping;            // IGA [no ref here]

bspline_mass_matrix.reinit (sparsity_bspline);
assemble_mass_matrix_diagonal(bspline_mass_matrix);
diagonal_of_mass_matrix.reinit (iga_handler.n_bspline);

// IGA MOD: Instead of extracting the diagonal of the mass matrix, we perform a mass lumping:
mass_lumping.reinit(iga_handler.n_bspline);
for (unsigned int i=0; i<mass_lumping.size(); ++i)
  mass_lumping[i] = 1;
bspline_mass_matrix.vmult(diagonal_of_mass_matrix, mass_lumping);

// IGA NEW: Boundary values
QGauss<dim-1>  boundary_quad(fe.degree+2);
std::map<types::global_dof_index,double> boundary_values;

typename std::map<types::boundary_id, const Function<dim>*> dirichlet_boundary;   // <- [ NEW ]
BoundaryValues<dim> boundary_funct;
dirichlet_boundary[0] = &boundary_funct;

iga_handler.project_boundary_values(dirichlet_boundary,
                                    boundary_quad,
                                    bspline_constraints);
```

---
### ↘️ make_grid()

FROM:
```cpp
 GridGenerator::hyper_cube(triangulation, -1, 1);
 triangulation.refine_global(7);
```
TO:
```cpp
 -nd-
```



