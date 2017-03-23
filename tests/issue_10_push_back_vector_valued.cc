#include <deal.II/lac/vector.h>
#include <deal.II/base/parsed_function.h>

#include <iostream>
#include <assert.h> 


/*!

@brief Test for issue https://github.com/alexanderzimmerman/nsb-pcm/issues/10

@detail 
    If the bug still exists, then this test will seg fault when calling parse_parameters.

*/
int main(int /*argc*/, char** /*argv*/)
{
    const unsigned int dim = 2;
    const unsigned int function_count = 4;
    const unsigned int vector_component_count = dim + 2;

    std::vector<std::shared_ptr<dealii::Functions::ParsedFunction<dim>>>
        function_pointers;

    dealii::ParameterHandler prm;
    dealii::Point<dim> point;
    dealii::Vector<double> vector_value(vector_component_count);

    for (unsigned int f = 0; f < function_count; ++f)
    {
        function_pointers.push_back(
            std::shared_ptr<dealii::Functions::ParsedFunction<dim>>(
                new dealii::Functions::ParsedFunction<dim>(
                    vector_component_count)));
            
        prm.enter_subsection("parsed_function_"+std::to_string(f));
        {
            function_pointers[f]->declare_parameters(prm, vector_component_count);

            function_pointers[f]->parse_parameters(prm);
        }
        prm.leave_subsection();

        function_pointers[f]->vector_value(point, vector_value);
            
            std::cout << "f(" << point << ") = " << vector_value << std::endl;
    }

    return 0;
}
