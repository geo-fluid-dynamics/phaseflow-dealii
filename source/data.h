#ifndef parameter_reader_h
#define parameter_reader_h

#include <deal.II/base/function.h>
#include <deal.II/base/parsed_function.h>

#include <iostream>
#include <fstream>
#include <assert.h>


/*!

@brief Contains data (parameters) class and related data structures.

@detail

    This will contain more parameters, but this file is a simplified test.

*/
namespace Data
{

    class Data
    {
        
        public:

            dealii::ParameterHandler prm;

            Data(dealii::ParameterHandler &_prm) : prm(_prm)
            {
                /*! Declare parameters. */
                this->declare();

                /*! Read parameters. */
                if (parameter_file != "")
                {
                    this->prm.read_input();    
                }

                this->get_data();
            }

            void write(const std::string parameter_file_path="used_parameters.prm") const;

        private:
            virtual void declare() const;
            virtual void get_data();
    };



    void Data::write(const std::string parameter_file_path)
    {
        // Print a log file of all the ParameterHandler parameters
        std::ofstream parameter_log_file(parameter_file_path);
        assert(parameter_log_file.good());
        this->prm.print_parameters(parameter_log_file, dealii::ParameterHandler::Text);
    }


    class MetaData : public Data
    {
        public:
            unsigned int dim;
            unsigned int boundary_count;
            unsigned int vector_component_count;

        private:
            void declare() const;
            void get_data();
    };

    void MetaData::declare() const
    {
        this->prm.enter_subsection("meta");
        
        this->prm.declare_entry(
            "dim", "2", dealii::Patterns::Integer(1, 3),
            "Set the number of dimensions in the spatial domain, i.e. the dim template parameter for most deal.II templates.");

        this->prm.declare_entry(
            "boundary_count", "2", dealii::Patterns::Integer(2, dealii::Patterns::max_int),
            "Set the number of boundaries on the domain."
            "This information is needed to declare the proper number of parsed boundary functions.");
            
        this->prm.declare_entry(
            "vector_component_count", "4",
            dealii::Patterns::Integer(),
            "Specify the number of components in the vector."
            "e.g. 4 for a 2D vector-valued velocity, scalar pressure, and scalar temperature.");

        this->prm.leave_subsection();
    }

    void MetaData::get_data(const dealii::ParameterHandler)
    {
        this->prm.enter_subsection("meta");

        this->dim = this->prm.get_integer("dim");  

        this->boundary_count = this->prm.get_integer("boundary_count");

        this->vector_component_count = this->prm.get_integer("vector_component_count");

        this->prm.leave_subsection();
    }


    /*!

    @brief This is a class for boundary conditions data.

    */
    template<int dim>
    class BoundaryConditionsData : public Data
    {
        public:
            std::vector<dealii::Functions::ParsedFunction<dim>> functions;

        private:
            virtual void declare() const;
            virtual void get_data();
    };

    template<int dim>
    void BoundaryConditionsData::declare() const
    {
        
        prm.enter_subsection("meta");

        unsigned int boundary_count = this->prm.get_integer("boundary_count");
        unsigned int vector_component_count = this->prm.get_integer("vector_component_count");

        prm.leave_subsection();


        prm.enter_subsection("boundary_conditions");

        /*! Declare parsed functions. */

        std::vector<dealii::Functions::ParsedFunction<dim>> functions(boundary_count);

        for (unsigned int b = 0; b < boundary_count; ++b)
        {
            prm.enter_subsection("parsed_function_"+std::to_string(b));
            {
                functions[b].declare_parameters(this->prm, vector_component_count);    
            }
            prm.leave_subsection();
        }

        prm.leave_subsection();

    }

    template<int dim>
    void BoundaryConditionsData::get_data()
    {

        this->prm.enter_subsection("meta");

        unsigned int boundary_count = this->prm.get_integer("boundary_count");

        this->prm.leave_subsection();

        
        this->functions.resize(boundary_count);


        this->prm.enter_subsection("boundary_conditions");

        for (unsigned int b = 0; b < boundary_count; ++b)
        {
            this->prm.enter_subsection("parsed_function_"+std::to_string(b));
            {
                this->boundary_functions[b].parse_parameters(this->prm);
            }
            this->prm.leave_subsection();
        }

        this->prm.leave_subsection();


        return data;
    }

    }
        

    template<int dim>
    class AllData
    {
        public:
            MetaData meta;
            BoundaryConditionsData<dim> boundary_conditions;

            AllData(dealii::ParameterHandler &_prm) 
                : prm(&_prm), meta(&_prm), boundary_conditions(&_prm)
            {}

            void write(const std::string parameter_file_path="used_parameters.prm") const;

        private:
            dealii::ParameterHandler prm;
            virtual void declare() const;
            virtual void get_data();
    };

    template<int dim>
    virtual void AllData::declare() const
    {}

    template<int dim>
    virtual void AllData::get_data() const
    {}

}

#endif