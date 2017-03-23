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
            unsigned int boundary_count;
            unsigned int vector_component_count;

            std::vector<std::shared_ptr<dealii::Functions::ParsedFunction<dim>>> function_pointers;

            void declare(dealii::ParameterHandler &prm) const;
            void get_data(dealii::ParameterHandler &prm);
    };

    template<int dim>
    void BoundaryConditionsData<dim>::declare(dealii::ParameterHandler &prm) const
    {
        prm.enter_subsection("boundary_conditions");

        /*! Declare parsed functions. */

        std::vector<dealii::Functions::ParsedFunction<dim>> 
            functions(this->boundary_count);

        for (unsigned int b = 0; b < boundary_count; ++b)
        {
            prm.enter_subsection("parsed_function_"+std::to_string(b));
            {
                functions[b].declare_parameters(prm, this->vector_component_count);    
            }
            prm.leave_subsection();
        }

        prm.leave_subsection();

    }

    template<int dim>
    void BoundaryConditionsData<dim>::get_data(dealii::ParameterHandler &prm)
    {
        
        prm.enter_subsection("boundary_conditions");

        for (unsigned int b = 0; b < this->boundary_count; ++b)
        {
            this->function_pointers.push_back(
                std::shared_ptr<dealii::Functions::ParsedFunction<dim>>(
                    new dealii::Functions::ParsedFunction<dim>(
                        this->vector_component_count)));

            prm.enter_subsection("parsed_function_"+std::to_string(b));
            {
                this->function_pointers[b]->parse_parameters(prm);
            }
            prm.leave_subsection();
        }

        prm.leave_subsection();

    }


    template<int dim>
    class AllData : public Data
    {
        public:
            unsigned int boundary_count;
            unsigned int vector_component_count;
            BoundaryConditionsData<dim> boundary_conditions;

            virtual void read(
                dealii::ParameterHandler &prm,
                const std::string parameter_file_path);

            void declare(dealii::ParameterHandler &prm) const;
            void get_data(dealii::ParameterHandler &prm);
    };

    template<int dim>
    void AllData<dim>::declare(dealii::ParameterHandler &prm) const
    {

        prm.declare_entry(
            "boundary_count", "4", dealii::Patterns::Integer(2, dealii::Patterns::Integer::max_int_value),
            "Set the number of boundaries on the domain."
            "This information is needed to declare the proper number of parsed boundary functions.");

        prm.declare_entry(
            "vector_component_count", "4",
            dealii::Patterns::Integer(),
            "Specify the number of components in the vector."
            "e.g. 4 for a 2D vector-valued velocity, scalar pressure, and scalar temperature.");

    }

    template<int dim>
    void AllData<dim>::get_data(dealii::ParameterHandler &prm)
    {
        this->boundary_count = prm.get_integer("boundary_count");

        this->vector_component_count = prm.get_integer("vector_component_count");
    }

    template<int dim>
    void AllData<dim>::read(
        dealii::ParameterHandler &prm,
        const std::string parameter_file_path)
    {
        /*! Declare parameters. */
        this->declare(prm);

        /*! Read parameters. */
        if (parameter_file_path != "")
        {
            prm.read_input(parameter_file_path);    
        }

        /*! Organize parameters in this class's data structure*/
        this->get_data(prm);

        /*! Now we have the required info to declare boundary conditions. */
        this->boundary_conditions.boundary_count = this->boundary_count;
        this->boundary_conditions.vector_component_count = this->vector_component_count;

        this->boundary_conditions.read(prm, parameter_file_path);

    }

}

#endif