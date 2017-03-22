#include <deal.II/lac/vector.h>

#include <iostream>
#include <assert.h> 

#include "peclet_data.h"

int main(int /*argc*/, char** /*argv*/)
{
    const unsigned int dim = 2;

    dealii::ParameterHandler prm;

    Peclet::AllData<dim> data;
    
    data.read(prm, "");
    
    dealii::Point<dim> point;
    dealii::Vector<double> vector_value(
        data.vector_component_count);

    std::cout << "Evaluating boundary functions at point " << point << std::endl;

    for (auto bfp : data.boundary_conditions.function_pointers)
    {
        bfp->vector_value(point, vector_value);
        std::cout << vector_value << std::endl;
    }
    
    return 0;
}