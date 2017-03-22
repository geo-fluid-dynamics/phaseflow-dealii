#ifndef _data_h
#define _data_h

#include <deal.II/base/parameter_handler.h>

#include <iostream>
#include <fstream>
#include <assert.h>


/*!

@brief This is a class for handling data/parameters.

@detail

    Parameter input file reading uses dealii::ParameterHandler.

@author Alexander G. Zimmerman <zimmerman@aices.rwth-aachen.de>

*/
class Data
{
    public:
        Data() {};
        void read(const std::string parameter_file_path="");

    private:
        virtual void declare(dealii::ParameterHandler &prm) const = 0;
        virtual void get_data(dealii::ParameterHandler &prm) = 0;
};

void Data::read(const std::string parameter_file_path)
{
    dealii::ParameterHandler prm;

    /*! Declare parameters. */
    this->declare(prm);

    /*! Read parameters. */
    if (parameter_file_path != "")
    {
        prm.read_input(parameter_file_path);    
    }

    this->get_data(prm);
}


/*!

@brief This derived class demonstrates how to use the Data class.

*/
class TestData : public Data
{
    public:
        bool pass;

    private:
        void declare(dealii::ParameterHandler &prm) const;
        void get_data(dealii::ParameterHandler &prm);
};

void TestData::declare(dealii::ParameterHandler &prm) const
{
    prm.enter_subsection("test_section");
    
    prm.declare_entry(
        "pass", "true", dealii::Patterns::Bool(),
        "This is for test/data.cc. Set to true to pass or false to fail.");

    prm.leave_subsection();
}

void TestData::get_data(dealii::ParameterHandler &prm)
{
    prm.enter_subsection("test_section");

    this->pass = prm.get_bool("pass");  

    prm.leave_subsection();
}

#endif