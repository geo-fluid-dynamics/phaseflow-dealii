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

    std::vector<std::shared_ptr<dealii::Functions::ParsedFunction<dim>>>
        function_pointers;

    dealii::ParameterHandler prm;

    for (unsigned int f = 0; f < function_count; ++f)
    {
        function_pointers.push_back(
            std::shared_ptr<dealii::Functions::ParsedFunction<dim>>(new dealii::Functions::ParsedFunction<dim>()));
            
        prm.enter_subsection("parsed_function_"+std::to_string(f));
        {
            function_pointers[f]->declare_parameters(prm);
        }
        prm.leave_subsection();
    }

    prm.read_input("default.prm", false, true);

    for (unsigned int f = 0; f < function_count; ++f)
    {
        prm.enter_subsection("parsed_function_"+std::to_string(f));
        {
            function_pointers[f]->parse_parameters(prm);
        }
        prm.leave_subsection();
    }

    /* Evaluate the functions at the default point. */
    dealii::Point<dim> point;
    double value;

    for (unsigned int f = 0; f < function_count; ++f)
    {
        function_pointers[f]->value(point, value);
        std::cout << "f(" << point << ") = " << value << std::endl;
    }

    return 0;
}
