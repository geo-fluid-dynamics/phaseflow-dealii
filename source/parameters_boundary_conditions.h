#ifndef parameters_boundary_conditions_h
#define parameters_boundary_conditions_h

#include <deal.II/base/function.h>
#include <deal.II/base/parsed_function.h>

#include <fstream>


/*!

@brief Contain methods and data structures related to boundary condition parameters.

@author Alexander G. Zimmerman <zimmerman@aices.rwth-aachen.de> 2017

*/
namespace Parameters
{

    /*!

    @brief Contain methods and data structures related to boundary condition parameters.

    @author Alexander G. Zimmerman <zimmerman@aices.rwth-aachen.de> 2017

    */
    namespace BoundaryConditions
    {

        
    




        /*!

        @brief Read input parameters for parsed boundary functions.

        @detail

            This will read more parameters, but this file is a simplified test.

        */
        template<int dim>
        BoundaryConditions read(
            const std::string parameter_file,
            const unsigned int vector_size,
            const unsigned int boundary_count)
        {
            BoundaryConditions bcs;
            
            read_parsed_functions(
                const std::string parameter_file,
                const unsigned int vector_size,
                const unsigned int boundary_count,
                std::vector<dealii::Functions::ParsedFunction<dim>> bcs.boundary_functions);

            return bcs;
        }

    }

}

#endif