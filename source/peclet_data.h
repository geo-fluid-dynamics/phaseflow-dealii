#ifndef _peclet_data_h
#define _peclet_data_h

#include "data.h"

#include <deal.II/base/function.h>
#include <deal.II/base/parsed_function.h>

namespace Peclet
{

    /*!

    @brief This is a class for boundary conditions data.

    */
    template<int dim>
    class BoundaryConditionsData : public Data
    {
        public:
            std::vector<dealii::Functions::ParsedFunction<dim>> functions;

        private:
            void declare(dealii::ParameterHandler &prm) const;
            void get_data(dealii::ParameterHandler &prm);
    };

    template<int dim>
    void BoundaryConditionsData<dim>::declare(dealii::ParameterHandler &prm) const
    {
        
        prm.enter_subsection("meta");

        unsigned int boundary_count = prm.get_integer("boundary_count");
        unsigned int vector_component_count = prm.get_integer("vector_component_count");

        prm.leave_subsection();


        prm.enter_subsection("boundary_conditions");

        /*! Declare parsed functions. */

        std::vector<dealii::Functions::ParsedFunction<dim>> functions(boundary_count);

        for (unsigned int b = 0; b < boundary_count; ++b)
        {
            prm.enter_subsection("parsed_function_"+std::to_string(b));
            {
                functions[b].declare_parameters(prm, vector_component_count);    
            }
            prm.leave_subsection();
        }

        prm.leave_subsection();

    }

    template<int dim>
    void BoundaryConditionsData<dim>::get_data(dealii::ParameterHandler &prm)
    {

        unsigned int boundary_count = prm.get_integer("boundary_count");
        
        this->functions.resize(boundary_count);

        prm.enter_subsection("boundary_conditions");

        for (unsigned int b = 0; b < boundary_count; ++b)
        {
            prm.enter_subsection("parsed_function_"+std::to_string(b));
            {
                this->boundary_functions[b].parse_parameters(prm);
            }
            prm.leave_subsection();
        }

        prm.leave_subsection();

    }


    template<int dim>
    class AllData
    {
        public:
            unsigned int boundary_count;
            unsigned int vector_component_count;
            BoundaryConditionsData<dim> boundary_conditions;

        private:
            void declare(dealii::ParameterHandler &prm) const;
            void get_data(dealii::ParameterHandler &prm);
    };

    template<int dim>
    void AllData<dim>::declare(dealii::ParameterHandler &prm) const
    {

        prm.declare_entry(
            "boundary_count", "2", dealii::Patterns::Integer(2, dealii::Patterns::Integer::max_int_value),
            "Set the number of boundaries on the domain."
            "This information is needed to declare the proper number of parsed boundary functions.");

        prm.declare_entry(
            "vector_component_count", "4",
            dealii::Patterns::Integer(),
            "Specify the number of components in the vector."
            "e.g. 4 for a 2D vector-valued velocity, scalar pressure, and scalar temperature.");

        this->boundary_conditions.declare(prm);

    }

    template<int dim>
    void AllData<dim>::get_data(dealii::ParameterHandler &prm)
    {
        this->boundary_count = prm.get_integer("boundary_count");

        this->vector_component_count = prm.get_integer("vector_component_count");

        this->boundary_conditions.get_data(prm);
    }

}

#endif