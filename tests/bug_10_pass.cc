#include <deal.II/base/parsed_function.h>

#include <iostream>
#include <assert.h> 


/*!

@brief Test for issue https://github.com/alexanderzimmerman/nsb-pcm/issues/10

*/
int main(int /*argc*/, char** /*argv*/)
{
    const unsigned int dim = 2;
    const unsigned int function_count = 4;

    std::vector<dealii::Functions::ParsedFunction<dim>> functions(function_count);

    dealii::ParameterHandler prm;

    for (unsigned int f = 0; f < function_count; ++f)
    {
        prm.enter_subsection("parsed_function_"+std::to_string(f));
        {
            functions[f].declare_parameters(prm);
        }
        prm.leave_subsection();
    }

    prm.read_input("default.prm", false, true);

    for (unsigned int f = 0; f < function_count; ++f)
    {
        prm.enter_subsection("parsed_function_"+std::to_string(f));
        {
            functions[f].parse_parameters(prm);
        }
        prm.leave_subsection();
    }


    /* Evaluate the functions at the default point. */
    dealii::Point<dim> point;
    double value;

    for (unsigned int f = 0; f < function_count; ++f)
    {
        functions[f].value(point, value);
        std::cout << "f(" << point << ") = " << value << std::endl;
    }

    return 0;
}
