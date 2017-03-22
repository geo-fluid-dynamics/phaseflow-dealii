#include <deal.II/lac/vector.h>

#include <iostream>
#include <assert.h> 

#include "data.h"

int main(int /*argc*/, char** /*argv*/)
{
    const unsigned int dim = 2;

    auto data = Data::read<dim>();
    
    dealii::Point<dim> point;
    dealii::Vector<double> vector_value(
        data.boundary_conditions.vector_component_count);

    for (auto bf : data.boundary_conditions.functions)
    {
        bf.vector_value(point, vector_value);
        std::cout << vector_value << std::endl;
    }
    
    return 0;
}