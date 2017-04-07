#include <iostream>

#include <deal.II/base/tensor.h>
#include <deal.II/base/numbers.h>

/*! This test confirms the equivalence of the tensor contractions v*T*w and w*T*v */
int main(int argc, char* argv[])
{
    double epsilon = 1.e-14;
   
    const unsigned int dim = 2;
   
    dealii::Tensor<1, dim> v({0.1419, 0.4218});
    
    dealii::Tensor<1, dim> w({0.9157, 0.7922});
    
    dealii::Tensor<2, dim> T;
    
    T[0][0] = 0.9595;
    T[0][1] = 0.0357;
    T[1][0] = 0.6557;
    T[1][1] = 0.8491;
 
    auto c = [](
        const dealii::Tensor<1, dim> _w,
        const dealii::Tensor<2, dim> _T,
        const dealii::Tensor<1, dim> _v)
    {
        return _v*_T*_w;
    };
    
    double vTw = c(v, T, w), wTv = c(w, T, v);
    
    assert(dealii::numbers::NumberTraits<double>::abs(vTw - wTv) > epsilon);
    
    return 0;
}

