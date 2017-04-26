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

#include "peclet_global_parameters.h"

namespace Peclet
{
  using namespace dealii;
    
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
    
    void setup_system();
    
    void assemble_system();
    
    void apply_boundary_values_and_constraints();
    
    void solve_linear_system();
    
    void step_newton();
    
    bool solve_nonlinear_problem();
    
    void set_time_step_size(double new_size);
    
    void step_time();

    void write_solution();

    Triangulation<dim> triangulation;

    FESystem<dim,dim> fe;
    
    const FEValuesExtractors::Vector velocity_extractor;
    
    const FEValuesExtractors::Scalar pressure_extractor;

    const FEValuesExtractors::Scalar temperature_extractor;
    
    DoFHandler<dim> dof_handler;

    ConstraintMatrix constraints;

    SparsityPattern sparsity_pattern;

    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    
    Vector<double> newton_residual;
    
    Vector<double> newton_solution;
    
    Vector<double> old_solution;
    
    Vector<double> old_newton_solution;

    Vector<double> system_rhs;
    
    double time;
    
    double time_step_size;
    
    unsigned int time_step_counter;

    unsigned int boundary_count;
    
    /*! These ID's label manifolds used for exact geometry */
    std::vector<unsigned int> manifold_ids;
        
    /*! These strings label types of manifolds used for exact geometry */
    std::vector<std::string> manifold_descriptors;
    
    Functions::ParsedFunction<dim> source_function;
    
    Functions::ParsedFunction<dim> initial_values_function;
    
    Functions::ParsedFunction<dim> boundary_function;
    
    Functions::ParsedFunction<dim> exact_solution_function;

    void append_verification_table();
    
    void write_verification_table();
    
    TableHandler verification_table;
    
    std::string verification_table_file_name = "verification_table.txt";
    
  };
  
  template<int dim>
  Peclet<dim>::Peclet()
    :
    fe(FE_Q<dim>(SCALAR_DEGREE + 1), dim, // velocity
       FE_Q<dim>(SCALAR_DEGREE), 1, // pressure
       FE_Q<dim>(SCALAR_DEGREE), 1), // temperature
    velocity_extractor(0),
    pressure_extractor(dim),
    temperature_extractor(dim + 1),
    dof_handler(this->triangulation),
    initial_values_function(dim + 1 + ENERGY_ENABLED),
    source_function(dim + 1 + ENERGY_ENABLED),
    exact_solution_function(dim + 1 + ENERGY_ENABLED)
  {}
  
  #include "peclet_grid.h"

  #include "peclet_system.h"

  #include "peclet_solve_nonlinear_problem.h"
  
  #include "peclet_step_time.h"
  
  #include "peclet_output.h"
  
  #include "peclet_verification.h"
  
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
    
    /* Attach manifolds for exact geometry 
    
    For now this only supports a single spherical manifold centered at the origin.
    
    */
    SphericalManifold<dim> spherical_manifold;
    
    for (unsigned int i = 0; i < manifold_ids.size(); i++)
    {
        if (manifold_descriptors[i] == "spherical")
        {
            this->triangulation.set_manifold(manifold_ids[i], spherical_manifold);      
        }
    }
    
    // Run initial refinement cycles
    
    this->triangulation.refine_global(this->params.refinement.initial_global_cycles);
    
    Refinement::refine_mesh_near_boundaries(
        this->triangulation,
        this->params.refinement.boundaries_to_refine,
        this->params.refinement.initial_boundary_cycles);
    
    // Initialize the linear system
    
    this->setup_system(); 

    this->time = 0.;
    
    this->set_time_step_size(this->params.time.initial_step_size);
    
    this->time_step_counter = 0;
    
    VectorTools::interpolate(
        this->dof_handler,
        this->initial_values_function,
        this->solution); 
    
    this->write_solution();
    
    for (; this->time_step_counter < this->params.time.max_steps; ++this->time_step_counter)
    { 
        if (this->time > (this->params.time.end - EPSILON))
        {
            break;
        }
        
        this->step_time();
        
        this->write_solution();

        if (this->params.verification.enabled)
        {
            this->append_verification_table();
        }
        
        if (this->params.time.stop_when_steady)
        {
            Vector<double> time_residual = this->solution; // There is evidently no Vector<Number> - Vector<Number> method.
            time_residual -= this->old_solution;
            
            double unsteadiness = time_residual.l2_norm()/this->solution.l2_norm();
            
            std::cout << "Unsteadiness, || w_{n+1} - w_n || / || w_{n+1} || = " << unsteadiness << std::endl;
            
            if (unsteadiness < this->params.time.steady_tolerance)
            {
                std::cout << "Reached steady state." << std::endl;
                break;
            }
            
        }
    
    } 
    
    /* Clean up. 
    
    Manifolds must be detached from Triangulations before leaving this scope.
    
    */
    this->triangulation.set_manifold(0);
    
  }
  
}
