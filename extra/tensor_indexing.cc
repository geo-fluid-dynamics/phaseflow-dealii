#include <iostream>

#include <deal.II/base/tensor.h>

int main(int argc, char* argv[])
{
    const unsigned int dim = 2;
   
    dealii::Tensor<2, dim> T;
    
    std::cout << "T = ";
    
    for (unsigned int i = 0; i < T.n_independent_components; ++i)
    {
        auto indices = T.unrolled_to_component_indices(i);
        
        std::cout << "T_{" << indices[0] + 1 << indices[1] + 1 << "}" << " ";
        
    }
    
    for (unsigned int i = 1; i < (dim + 1); ++i)
    {
        for (unsigned int j = 1; j < (dim + 1); ++j)
        {
            T[i - 1][j - 1] = 10*i + j;
        }
    }

    std::cout << std::endl << std::endl << "T = " << T << std::endl;
    
    return 0;
}
