#ifndef peclet_parameters_h
#define peclet_parameters_h

#include <vector>
#include <iostream>
#include <fstream>
#include <functional>

#include "my_parameter_handler.h"

#include "peclet_global_parameters.h"

/*
    
    @brief Encapsulates parameter handling and paramter input file handling.

    @detail
    
        Originally the ParameterReader from deal.II's step-26 was used;
        but that was too prohibitive. The approach in this file isolates the details
        of the input file handling from the rest of the program,
        allowing for rich data structures and simplifying the user code.

        The goal is to allow one to run the program with different parameters,
        without having to recomplile the program.

        It is valuable to structure the parameter data as done here, to greatly
        simplify the writing and debugging of the code.

        This also makes it simple to insantiate a PDE Model in a user program 
        and then change it's parameters directly without having to use any intermediate text files.

        Every parameter name appears three times in this file. Once for the data structure,
        once for declaring to the ParameterHandler class, and once for parsing the input file.
    
    @author A. Zimmerman <zimmerman@aices.rwth-aachen.de>
    
*/

namespace Peclet
{
    namespace Parameters
    {   

        using namespace dealii;

        struct Meta
        {
            unsigned int dim;
        };

        struct PhysicalModel
        {
            std::vector<double> gravity;
        };
        
        struct Geometry
        {
            unsigned int dim;
            std::string grid_name;
            std::vector<double> sizes;
            std::vector<double> transformations;
        };
        
        struct BoundaryConditions
        {
            std::vector<std::string> strong_mask;
        };
        
        struct AdaptiveRefinement
        {
            unsigned int initial_cycles;
            unsigned int max_level;
            unsigned int max_cells;
            unsigned int interval;
            unsigned int cycles_at_interval;
            double refine_fraction;
            double coarsen_fraction;
        };
            
        struct Refinement
        {
            unsigned int initial_global_cycles;
            unsigned int initial_boundary_cycles;
            std::vector<unsigned int> boundaries_to_refine;
            AdaptiveRefinement adaptive;
        };
        
        struct Time
        {
            double end;
            double initial_step_size;
            double min_step_size;
            double max_step_size;
            unsigned int max_steps;
            bool stop_when_steady;
            double steady_tolerance;
        };

        struct IterativeSolver
        {
            std::string method;
            unsigned int max_iterations;
            double tolerance;
        };
        
        struct Output
        {
            bool write_solution_vtk;
        };
        
        struct Verification
        {
            bool enabled;
        };
        
        struct StructuredParameters
        {
            Meta meta;
            PhysicalModel physics;
            BoundaryConditions boundary_conditions;
            Geometry geometry;
            Refinement refinement;
            Time time;
            IterativeSolver nonlinear_solver;
            Output output;
            Verification verification;
        };    

        template<int dim>
        void declare(ParameterHandler &prm)
        {
            
            prm.enter_subsection("meta");
            {
                prm.declare_entry("dim", std::to_string(dim), Patterns::Integer(1, 3));
            }
            prm.leave_subsection();

            
            prm.enter_subsection("physics");
            {
                prm.declare_entry("gravity", "0., -1, 0.", Patterns::List(Patterns::Double()));
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("source_function");
            {
                Functions::ParsedFunction<dim>::declare_parameters(prm, dim + 2);    
            }
            prm.leave_subsection();  

            
            prm.enter_subsection ("geometry");
            {
                    
                prm.declare_entry("grid_name", "hyper_rectangle",
                     Patterns::Selection("hyper_rectangle | hyper_shell"),
                     "Select the name of the geometry and grid to generate.");
                     
                prm.declare_entry("sizes", "0., 0., 1., 1.",
                    Patterns::List(Patterns::Double(0.)),
                    "Set the sizes for the grid's geometry.");
                              
            }
            prm.leave_subsection ();
            

            prm.enter_subsection ("initial_values");
            {
                Functions::ParsedFunction<dim>::declare_parameters(prm, dim + 2); 
            }
            prm.leave_subsection ();
            
            
            prm.enter_subsection ("boundary_conditions");
            {
                /* It was originally attempted to make this parameter a list of lists,
                but Patterns::List does not appear to allow for that.
                So instead we have one string that we'll parse manually.
                */
                prm.declare_entry(
                    "strong_mask",
                    "velocity",
                    Patterns::List(Patterns::Selection("velocity | pressure | temperature")),
                    "The parsed functions will only be applied as strong boundary conditions to components included in the mask.");

            }
            prm.leave_subsection ();
            
            
            prm.enter_subsection ("boundary_conditions");
            {
                /* It was originally attempted to make this parameter a list of lists,
                but Patterns::List does not appear to allow for that.
                So instead we have one string that we'll parse manually.
                */
                prm.declare_entry(
                    "strong_mask",
                    "velocity",
                    Patterns::List(Patterns::Selection("velocity | pressure | temperature")),
                    "The parsed functions will only be applied as strong boundary conditions to components included in the mask.");

            }
            prm.leave_subsection ();
            
            
            prm.enter_subsection ("refinement");
            {
                prm.declare_entry("initial_global_cycles", "2",
                    Patterns::Integer(),
                    "Initially globally refine the grid this many times "
                    "without using any error measure");
                    
                prm.declare_entry("initial_boundary_cycles", "0",
                    Patterns::Integer(),
                    "Initially refine the grid this many times"
                    "near the boundaries that are listed for refinement");
                    
                prm.declare_entry("boundaries_to_refine", "0",
                    Patterns::List(Patterns::Integer()),
                    "Refine cells that contain these boundaries");
                
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection ("time");
            {
                prm.declare_entry("end", "0.",
                    Patterns::Double(0.),
                    "End the time-dependent simulation once this time is reached.");
                
                prm.declare_entry("initial_step_size", "1.",
                    Patterns::Double(0.),
                    "Begin with this time step size.");
                
                prm.declare_entry("min_step_size", "1.",
                    Patterns::Double(0.),
                    "Minimum step size for adaptive time steppinig.");
                    
                prm.declare_entry("max_step_size", "1.",
                    Patterns::Double(0.),
                    "Maximum step size for adaptive time steppinig.");
                    
                prm.declare_entry("max_steps", "1000000",
                    Patterns::Integer(0),
                    "Maximum number of time steps.");
                    
                prm.declare_entry("stop_when_steady", "false", Patterns::Bool());
                
                prm.declare_entry("steady_tolerance", "1.e-8", Patterns::Double(0.));
                    
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("nonlinear_solver");
            {
                prm.declare_entry("method", "Newton",
                     Patterns::Selection("Newton"));
                     
                prm.declare_entry("max_iterations", "50",
                    Patterns::Integer(0));
                    
                prm.declare_entry("tolerance", "1e-9",
                    Patterns::Double(0.));
                    
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("output");
            {
                prm.declare_entry("write_solution_vtk", "true", Patterns::Bool());
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("verification");
            {
                prm.declare_entry("enabled", "false", Patterns::Bool());

                prm.enter_subsection("exact_solution_function");
                {
                    Functions::ParsedFunction<dim>::declare_parameters(prm, dim + 2);    
                }
                prm.leave_subsection();
            }
            prm.leave_subsection();

        }   

        
        Meta read_meta_parameters(const std::string parameter_file="")
        {
            Meta mp;
            
            ParameterHandler prm;
            declare<2>(prm);
            
            if (parameter_file != "")
            {
                prm.parse_input(parameter_file);    
            }
            
            prm.enter_subsection("meta");
            {
                mp.dim = prm.get_integer("dim");  
            }
            prm.leave_subsection();

            return mp;
        }
        
        template <int dim>
        StructuredParameters read(
                const std::string parameter_file,
                Functions::ParsedFunction<dim> &source_function,
                Functions::ParsedFunction<dim> &initial_values_function,
                Functions::ParsedFunction<dim> &exact_solution_function)
        {

            StructuredParameters params;
            
            ParameterHandler prm;
            Parameters::declare<dim>(prm);

            if (parameter_file != "")
            {
                prm.parse_input(parameter_file);
            }
            
            // Print a log file of all the ParameterHandler parameters
            std::ofstream parameter_log_file("used_parameters.prm");
            assert(parameter_log_file.good());
            prm.print_parameters(parameter_log_file, ParameterHandler::Text);
            
            prm.enter_subsection("physics");
            {
                params.physics.gravity = MyParameterHandler::get_vector<double>(prm, "gravity");
            }
            prm.leave_subsection();
            
            prm.enter_subsection("geometry");
            {
                params.geometry.grid_name = prm.get("grid_name");
                params.geometry.sizes = MyParameterHandler::get_vector<double>(prm, "sizes");
            }
            prm.leave_subsection();

            
            prm.enter_subsection("source_function");
            {
                source_function.parse_parameters(prm);
            }
            prm.leave_subsection();
                
            
            prm.enter_subsection("verification");
            {
                
                params.verification.enabled = prm.get_bool("enabled");
                        
                prm.enter_subsection("exact_solution_function");
                {
                    exact_solution_function.parse_parameters(prm);    
                }
                prm.leave_subsection();
                
            }
            prm.leave_subsection();

            
            prm.enter_subsection("initial_values");
            {               

                 initial_values_function.parse_parameters(prm);
                 
            }    
            prm.leave_subsection();
            
            
            prm.enter_subsection ("boundary_conditions");
            {
                params.boundary_conditions.strong_mask = 
                    MyParameterHandler::get_vector<std::string>(prm, "strong_mask");
            }
            prm.leave_subsection ();
            
            
            prm.enter_subsection("refinement");
            {
                
                params.refinement.initial_global_cycles = prm.get_integer("initial_global_cycles");
                params.refinement.initial_boundary_cycles = prm.get_integer("initial_boundary_cycles");
                params.refinement.boundaries_to_refine = 
                    MyParameterHandler::get_vector<unsigned int>(prm, "boundaries_to_refine");
                
            }
            prm.leave_subsection();
            
            
            prm.enter_subsection("time");
            {
                params.time.end = prm.get_double("end");
                params.time.initial_step_size = prm.get_double("initial_step_size");
                params.time.min_step_size = prm.get_double("min_step_size");
                params.time.max_step_size = prm.get_double("max_step_size");
                params.time.max_steps = prm.get_integer("max_steps");
                params.time.stop_when_steady = prm.get_bool("stop_when_steady");
                params.time.steady_tolerance = prm.get_double("steady_tolerance");
            }    
            prm.leave_subsection();
            
            
            prm.enter_subsection("nonlinear_solver");
            {
                params.nonlinear_solver.method = prm.get("method");
                params.nonlinear_solver.max_iterations = prm.get_integer("max_iterations");
                params.nonlinear_solver.tolerance = prm.get_double("tolerance");
            }    
            prm.leave_subsection(); 
            
            
            prm.enter_subsection("output");
            {
                params.output.write_solution_vtk = prm.get_bool("write_solution_vtk");
            }
            prm.leave_subsection();
            
            return params;
        }

    }    

}

#endif
