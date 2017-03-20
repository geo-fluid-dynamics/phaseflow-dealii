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

#include "extrapolated_field.h"
#include "my_grid_generator.h"
#include "fe_field_tools.h"
#include "output.h"

#include "peclet_parameters.h"

namespace Peclet
{
  using namespace dealii;

  const unsigned int SCALAR_DEGREE = 1; /*! @todo: Expose to ParameterHandler */
  const unsigned int VECTOR_DEGREE = 2; /*! @todo: Expose to ParameterHandler */
  const unsigned int MAX_NEWTON_ITERATIONS = 100; /*! @todo: Expose to ParameterHandler */
  const double NEWTON_TOLERANCE = 1.e-10; /*! @todo: Expose to ParameterHandler */
  
  const double RAYLEIGH_NUMBER = 1.e6; /*! @todo: Expose to ParameterHandler */
  const double PRANDTL_NUMBER = 0.71; /*! @todo: Expose to ParameterHandler */
  const double REYNOLDS_NUMBER = 1.;

  const double LIQUID_DYNAMIC_VISOCITY = 1.; /*! @todo: Expose to ParameterHandler */
  const double SOLID_CONDUCTIVITY = 1.; /*! @todo: Expose to ParameterHandler */
  const double LIQUID_CONDUCTIVITY = 2.; /*! @todo: Expose to ParameterHandler */
  
  const Tensor<1, 3> GRAVITY({0., -1., 0.}); /*! @todo: Expose to ParameterHandler */

  const double EPSILON = 1.e-14;
  
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
    void adaptive_refine();
    void setup_system(bool quiet = false);
    void assemble_system();
    void apply_boundary_values_and_constraints();
    void step_time(bool quiet = false);

    SolverStatus solve_linear_system(bool quiet = false);

    void write_solution();

    Triangulation<dim> triangulation;

    FESystem<dim,dim> fe;
    
    DoFHandler<dim> dof_handler;

    ConstraintMatrix constraints;

    SparsityPattern sparsity_pattern;

    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       old_solution;
    Vector<double>       newton_solution;
    Vector<double>       old_newton_solution;

    Vector<double>       system_rhs;

    double               time;
    double               time_step_size;
    unsigned int         time_step_counter;
    bool final_time_step;
    bool output_this_step;
    SolverStatus solver_status;

    Point<dim> spherical_manifold_center;
    
    std::vector<unsigned int> manifold_ids;
    std::vector<std::string> manifold_descriptors;
    
    Function<dim>* velocity_function;
    Function<dim>* source_function;
    std::vector<Function<dim>*> boundary_functions;
    Function<dim>* initial_values_function;
    Function<dim>* exact_solution_function;
    
    void append_verification_table();
    void write_verification_table();
    TableHandler verification_table;
    std::string verification_table_file_name = "verification_table.txt";
    
  };
  
  template<int dim>
  Peclet<dim>::Peclet()
    :
    fe(FE_Q<dim>(VECTOR_DEGREE), dim, // velocity
       FE_Q<dim>(SCALAR_DEGREE), 2),  // pressure and temperature
    dof_handler(this->triangulation)
  {}
  
  #include "peclet_grid.h"

  #include "peclet_system.h"

  #include "peclet_step_time.h"
  
  template<int dim>
  void Peclet<dim>::write_solution()
  {
      
    if (this->params.output.write_solution_vtk)
    {
        Output::write_solution_to_vtk(
            "solution-"+Utilities::int_to_string(this->time_step_counter)+".vtk",
            this->dof_handler,
            this->solution);    
    }
    
  }
  
  template<int dim>
  void Peclet<dim>::append_verification_table()
  {
    assert(this->params.verification.enabled);
    
    this->exact_solution_function->set_time(this->time);
    
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    
    VectorTools::integrate_difference(
        this->dof_handler,
        this->solution,
        *this->exact_solution_function,
        difference_per_cell,
        QGauss<dim>(dim + 1),
        VectorTools::L2_norm);
        
    double L2_norm_error = difference_per_cell.l2_norm();
    
    VectorTools::integrate_difference(
        this->dof_handler,
        this->solution,
        *this->exact_solution_function,
        difference_per_cell,
        QGauss<dim>(dim + 1),
        VectorTools::L1_norm);
        
    double L1_norm_error = difference_per_cell.l1_norm();
        
    this->verification_table.add_value("time_step_size", this->time_step_size);
    this->verification_table.add_value("time", this->time);
    this->verification_table.add_value("cells", this->triangulation.n_active_cells());
    this->verification_table.add_value("dofs", this->dof_handler.n_dofs());
    this->verification_table.add_value("L1_norm_error", L1_norm_error);
    this->verification_table.add_value("L2_norm_error", L2_norm_error);
    
  }
  
  template<int dim>
  void Peclet<dim>::write_verification_table()
  {
    const int precision = 14;
    
    this->verification_table.set_precision("time", precision);
    this->verification_table.set_scientific("time", true);
    
    this->verification_table.set_precision("time_step_size", precision);
    this->verification_table.set_scientific("time_step_size", true);
    
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
    Functions::ParsedFunction<dim> parsed_velocity_function(dim),
        parsed_source_function(dim + 2),
        parsed_boundary_function(dim + 2),
        parsed_initial_values_function(dim + 2),
        parsed_exact_solution_function(dim + 2); 
    
    this->params = Parameters::read<dim>(
        parameter_file,
        parsed_source_function,
        parsed_boundary_function,
        parsed_exact_solution_function,
        parsed_initial_values_function);
    
    this->create_coarse_grid();
    
    this->source_function = &parsed_source_function;
    this->exact_solution_function = &parsed_exact_solution_function;
    
    /*
    Generalizing the handling of auxiliary functions is complicated. In most cases one should be able to use a ParsedFunction, but the generality of Function<dim>* allows for a standard way to account for any possible derived class of Function<dim>. 
    For example this allows for....
        - an optional initial values function that interpolates an old solution loaded from disk. 
        - flexibily implementing general boundary conditions
    */

    // Initial values function
    
    Triangulation<dim> field_grid;
    DoFHandler<dim> field_dof_handler(field_grid);
    Vector<double> field_solution;
    
    if (this->params.initial_values.function_name != "interpolate_old_field")
    { // This will write files that need to exist.
        this->setup_system(true);
        FEFieldTools::save_field_parts(this->triangulation, this->dof_handler, this->solution); 
    }
    
    FEFieldTools::load_field_parts(
        field_grid,
        field_dof_handler,
        field_solution,
        this->fe);
    
    MyFunctions::ExtrapolatedField<dim> field_function(
        field_dof_handler,
        field_solution);
    

    if (this->params.initial_values.function_name == "interpolate_old_field")
    {
        this->initial_values_function = &field_function;                      
    }
    else if (this->params.initial_values.function_name == "parsed")
    { 
        this->initial_values_function = &parsed_initial_values_function;
    }
    
    // Boundary condition functions
    
    unsigned int boundary_count = this->params.boundary_conditions.implementation_types.size();
    
    assert(params.boundary_conditions.function_names.size() == boundary_count);

    std::vector<ConstantFunction<dim>> constant_functions;
    
    for (unsigned int boundary = 0; boundary < boundary_count; boundary++)
    {
        std::string boundary_type = this->params.boundary_conditions.implementation_types[boundary];
        std::string function_name = this->params.boundary_conditions.function_names[boundary];
        
        if (function_name == "constant")
        {
            double value = this->params.boundary_conditions.function_double_arguments.front();
            this->params.boundary_conditions.function_double_arguments.pop_front();
            constant_functions.push_back(ConstantFunction<dim>(value));
        }
    }

    /*! Organize boundary functions to simplify application during the time loop */
    
    unsigned int constant_function_index = 0;
    
    for (unsigned int boundary = 0; boundary < boundary_count; boundary++)        
    {
        std::string boundary_type = this->params.boundary_conditions.implementation_types[boundary];
        std::string function_name = this->params.boundary_conditions.function_names[boundary];

        if (function_name == "constant")
        {
            assert(constant_function_index < constant_functions.size());
            this->boundary_functions.push_back(&constant_functions[constant_function_index]);
            constant_function_index++;
        }
        else if (function_name == "parsed")
        {
            this->boundary_functions.push_back(&parsed_boundary_function);
        }
        
    }

    // Attach manifolds
    
    assert(dim < 3); // @todo: 3D extension: For now the CylindricalManifold is being ommitted.
        // deal.II makes is impractical for a CylindricalManifold to exist in 2D.
    SphericalManifold<dim> spherical_manifold(this->spherical_manifold_center);
    
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

    Vector<double> tmp;
    Vector<double> forcing_terms;
    
    // Iterate
    unsigned int pre_refinement_step = 0;
    
start_time_iteration:

    tmp.reinit(this->solution.size());

    VectorTools::interpolate(this->dof_handler,
                             *this->initial_values_function,
                             this->old_solution); 
    
    this->solution = this->old_solution;
    this->time_step_counter = 0;
    this->time = 0;
    
    this->time_step_size = this->params.time.step_size;
    if (this->time_step_size < EPSILON)
    {
        this->time_step_size = this->params.time.end_time/
            pow(2., this->params.time.global_refinement_levels);
    }
    double Delta_t = this->time_step_size;
    
    this->write_solution();

    this->final_time_step = false;
    this->output_this_step = true;

    double epsilon = 1e-14;
    
    do
    {
        ++this->time_step_counter;
        this->time = Delta_t*time_step_counter; // Incrementing the time directly would accumulate errors
        
        // Set some flags that will control output for this step.
        this->final_time_step = this->time > this->params.time.end_time - epsilon;
        
        bool at_interval = false;
        
        if (this->params.output.time_step_interval == 1)
        {
            at_interval = true;
        }
        else if (this->params.output.time_step_interval != 0)
        {
            if ((this->time_step_counter % this->params.output.time_step_interval) == 0)
            {
                at_interval = true;
            }
        }
        else
        {
            at_interval = false;
        }
        
        if (at_interval)
        {
            this->output_this_step = true;
        }
        else
        {
            this->output_this_step = false;
        }
        
        /*!
         Run the time step
         */

        this->step_time();
        
        if ((this->time_step_counter == 1) &&
            (pre_refinement_step < this->params.refinement.adaptive.initial_cycles))
        {
            this->adaptive_refine();
            ++pre_refinement_step;

            tmp.reinit(this->solution.size());

            std::cout << std::endl;

            goto start_time_iteration;
        }
        else if ((this->time_step_counter > 0) 
                 && (params.refinement.adaptive.interval > 0)  // % 0 (mod 0) is undefined
                 && (this->time_step_counter % params.refinement.adaptive.interval == 0))
        {
            for (unsigned int cycle = 0;
                 cycle < params.refinement.adaptive.cycles_at_interval; cycle++)
            {
                this->adaptive_refine();
            }
            tmp.reinit(this->solution.size());
            
        }
        
        this->old_solution = this->solution;
        
    } while (!this->final_time_step);
    
    // Save data for FEFieldFunction so that it can be loaded for initialization
    FEFieldTools::save_field_parts(triangulation, dof_handler, solution);
    
    // Write error table
    if (this->params.verification.enabled)
    {
        this->write_verification_table();
    }

    // Cleanup
    this->triangulation.set_manifold(0);
    
  }
}
