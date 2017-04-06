 /*
 * @brief nsb-pcm peclet solves the Navier-Stokes-Boussinesq equations coupled with phase-change
 *
 * @detail
 *
 *  Based on peclet, which was based on deal.II Tutorial 26 by Wolfgang Bangerth, Texas A&M University, 2013
 *
 *  Some of the more notable extensions include:
 *  - Step time with Newton iterations; each Newton iteration requires a linear system solve
 *  - Assembles Newton linearized differential of the Navier-Stokes-Boussinesq equation with latent heat of phase-change
 *  - Supports time-dependent non-zero Dirichlet and Neumann boundary condition
 *  - Re-designed parmameter handling
 *  - Generalized boundary condition handling via the parameter input file
 *  - Writes FEFieldFunction to disk, and can read it from disk to initialize a restart
 *  - Extended the FEFieldFunction class for extrapolation
 *  - Added verification via Method of Manufactured Solutions (MMS) with error table based on approach from Tutorial 7
 *  - Added test suite using ctest and the standard deal.II approach
 *  - Added a parameteric sphere-cylinder grid
 *  - Added a boundary grid refinement routine
 *
 * @author Alexander Zimmerman <zimmerman@aices.rwth-aachen.de>, RWTH AAchen University, 2016
 */
#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/table_handler.h>

#include <iostream>
#include <functional>

#include <assert.h> 
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/lac/sparse_direct.h>

#include "extrapolated_field.h"
#include "my_grid_generator.h"
#include "fe_field_tools.h"
#include "output.h"

#include "peclet_parameters.h"

#include "global_parameters.h"

namespace Peclet
{
  using namespace dealii;
  
  struct SolverStatus
  {
      unsigned int last_step;
  };
  
  template<int dim>
  class Peclet
  {
  public:
  
    Peclet();
    Parameters::StructuredParameters params;
    void init(std::string file_path);
    void run(const std::string parameter_file = "");

  private:

    void create_coarse_grid();
    void setup_system(bool quiet = false);
    void assemble_system();
    void apply_boundary_values_and_constraints();
    void step_newton();
    void solve_nonlinear_problem(bool quiet = false);

    SolverStatus solve_linear_system(bool quiet = false);

    void write_solution();

    Triangulation<dim> triangulation;

    FESystem<dim,dim> fe;
    
    DoFHandler<dim> dof_handler;

    ConstraintMatrix constraints;

    SparsityPattern sparsity_pattern;

    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       newton_residual;
    Vector<double>       old_newton_solution;

    Vector<double>       system_rhs;

    SolverStatus solver_status;

    unsigned int boundary_count;
    
    Functions::ParsedFunction<dim> initial_values_function;
    Functions::ParsedFunction<dim> source_function;
    Functions::ParsedFunction<dim> exact_solution_function;

    void append_verification_table();
    void write_verification_table();
    TableHandler verification_table;
    std::string verification_table_file_name = "verification_table.txt";
    
  };
  
  template<int dim>
  Peclet<dim>::Peclet()
    :
    fe(FE_Q<dim>(VECTOR_DEGREE), dim, // velocity
       FE_Q<dim>(SCALAR_DEGREE), 1),  // pressure
    dof_handler(this->triangulation),
    initial_values_function(dim + 1),
    source_function(dim + 1),
    exact_solution_function(dim + 1)
  {}
  
  #include "peclet_grid.h"

  #include "peclet_system.h"

  #include "peclet_solve_nonlinear_problem.h"
  
  template<int dim>
  void Peclet<dim>::write_solution()
  {
      
    if (this->params.output.write_solution_vtk)
    {
        Output::write_solution_to_vtk(
            "solution.vtk",
            this->dof_handler,
            this->solution);    
    }
    
  }
  
  template<int dim>
  void Peclet<dim>::append_verification_table()
  {
    assert(this->params.verification.enabled);
    
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    
    VectorTools::integrate_difference(
        this->dof_handler,
        this->solution,
        this->exact_solution_function,
        difference_per_cell,
        QGauss<dim>(dim + 1),
        VectorTools::L2_norm);
        
    double L2_norm_error = difference_per_cell.l2_norm();
    
    VectorTools::integrate_difference(
        this->dof_handler,
        this->solution,
        this->exact_solution_function,
        difference_per_cell,
        QGauss<dim>(dim + 1),
        VectorTools::L1_norm);
        
    double L1_norm_error = difference_per_cell.l1_norm();

    this->verification_table.add_value("cells", this->triangulation.n_active_cells());
    this->verification_table.add_value("dofs", this->dof_handler.n_dofs());
    this->verification_table.add_value("L1_norm_error", L1_norm_error);
    this->verification_table.add_value("L2_norm_error", L2_norm_error);
    
  }
  
  template<int dim>
  void Peclet<dim>::write_verification_table()
  {
    const int precision = 14;

    this->verification_table.set_precision("cells", precision);
    this->verification_table.set_scientific("cells", true);
    
    this->verification_table.set_precision("dofs", precision);
    this->verification_table.set_scientific("dofs", true);
    
    this->verification_table.set_precision("L2_norm_error", precision);
    this->verification_table.set_scientific("L2_norm_error", true);
    
    this->verification_table.set_precision("L1_norm_error", precision);
    this->verification_table.set_scientific("L1_norm_error", true);
    
    std::ofstream out_file(this->verification_table_file_name, std::fstream::app);
    assert(out_file.good());
    this->verification_table.write_text(out_file);
    out_file.close(); 
  }
  
  template<int dim>
  void Peclet<dim>::run(const std::string parameter_file)
  {
    
    // Clean up the files in the working directory

    if (this->params.verification.enabled)
    {
        std::remove(this->verification_table_file_name.c_str());
    }
    
    /*
    Working with deal.II's Function class has been interesting, and I'm 
    sure many of my choices are unorthodox. The most important lesson learned has been that 
    a Function<dim>* can point to any class derived from Function<dim>. The general design pattern
    then is to instantitate all of the functions that might be needed, and then to point to the ones
    actually being used.
    */
    
    this->params = Parameters::read<dim>(
        parameter_file,
        this->source_function,
        this->initial_values_function,
        this->exact_solution_function);
    
    this->create_coarse_grid();
    
    // Run initial refinement cycles
    
    this->triangulation.refine_global(this->params.refinement.initial_global_cycles);
    
    // Initialize the linear system
    
    this->setup_system(); 

    VectorTools::interpolate(this->dof_handler,
                             this->initial_values_function,
                             this->solution); 
    
    this->solve_nonlinear_problem();
    
    this->write_solution();
    
  }
  
}
