#include <assert.h> 

#include "data.h"

int main(int /*argc*/, char** /*argv*/)
{
    
    dealii::ParameterHandler prm;

    TestData data;
    
    data.read(prm, "");

    assert(data.pass);

    return 0;
}