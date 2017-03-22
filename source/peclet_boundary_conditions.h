#ifndef peclet_boundary_conditions_h
#define peclet_boundary_conditions_h


/*!
@brief Read input parameters for parsed boundary functions.

@detail 
    
    This is separate from the reading of other parameters, because the number of boundary functions is not known until after generating the coarse grid.

@author Alexander G. Zimmerman <zimmerman@aices.rwth-aachen.de> 2017

*/
template<int dim>
void Peclet<dim>::read_parsed_boundary_function_inputs(
    const std::string parameter_file,
    std::vector<Functions::ParsedFunction<dim>> &boundary_functions)
{
    
    ParameterHandler prm;

    prm.enter_subsection("boundary_conditions");

    for (unsigned int b = 0; b < BOUNDARY_COUNT; ++b)
    {
        prm.enter_subsection("parsed_function_"+std::to_string(b));
        {
            boundary_functions[b].declare_parameters(prm, dim + 2);    
        }
        prm.leave_subsection();
    }

    prm.leave_subsection();
        

    if (parameter_file != "")
    {
        prm.read_input(parameter_file);    
    }
    
    // Print a log file of all the ParameterHandler parameters
    std::ofstream parameter_log_file("used_boundary_parameters.prm");
    assert(parameter_log_file.good());
    prm.print_parameters(parameter_log_file, ParameterHandler::Text);
    
    prm.enter_subsection("boundary_conditions");
    {

        for (unsigned int b = 0; b < BOUNDARY_COUNT; ++b)
        {
            prm.enter_subsection("parsed_function_"+std::to_string(b));
            {
                boundary_functions[b].parse_parameters(prm);
            }
            prm.leave_subsection();
        }
    }    
    prm.leave_subsection();

}

#endif