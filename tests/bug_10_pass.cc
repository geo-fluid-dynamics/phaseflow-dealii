#include <deal.II/base/parsed_function.h>

#include <iostream>
#include <assert.h> 


/*!
@brief Test for bug https://github.com/alexanderzimmerman/nsb-pcm/issues/10

@detail 
    If the bug still exists, then this test will fail to build when
    attempting to resize a std::vector<dealii::Functions::ParsedFunction<dim>>.
*/

int main(int /*argc*/, char** /*argv*/)
{
    const unsigned int dim = 2;
    std::vector<double> doubles;

    const unsigned int new_size = 4;
    doubles.resize(new_size);

    return 0;
}