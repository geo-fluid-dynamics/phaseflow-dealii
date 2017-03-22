#include <assert.h> 

#include "data.h"

int main(int /*argc*/, char** /*argv*/)
{
    TestData data;

    data.read();

    assert(data.pass);

    return 0;
}