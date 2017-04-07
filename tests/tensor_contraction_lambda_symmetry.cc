#include <iostream>

#include <deal.II/base/tensor.h>
#include <deal.II/base/numbers.h>

/*! This test confirms the equivalence of the tensor contractions v*T*w and w*T*v */
int main(int argc, char* argv[])
{
    const unsigned int dim = 2;
    
    double epsilon = 1.e-14;
   
    dealii::Tensor<1, dim> v;
    dealii::Tensor<1, dim> w;
    dealii::Tensor<2, dim> T;
    
    for (unsigned int i = 0; i < dim; ++i)
    {
        v[i] = i + 1;
        
        w[i] = 2.*(i + 1);
        
        for (unsigned int j = 0; j < dim; ++j)
        {
            T[i][j] = 10*(i + 1) + j + 1;
        }
    }
    
    auto c = [](
        const Tensor<1, dim> _w,
        const Tensor<2, dim> _T,
        const Tensor<1, dim> _v)
    {
        return _w*_T*_v;
    };
    
    double vTw = c(v, T, w), wTv = c(w, T, v);
    
    assert(dealii::numbers::NumberTraits<double>::abs(vTw - wTv) < epsilon);
    
    return 0;
}

