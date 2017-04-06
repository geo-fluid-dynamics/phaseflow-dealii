#ifndef global_parameters_h
#define global_parameters_h

/*! Set global parameters.

Most of these aren't actually required at compile time, and should ideally be exposed to ParameterHandler.
This is just a temporary solution.

*/
const unsigned int SCALAR_DEGREE = 1; /*! @todo: Expose to ParameterHandler */

const unsigned int VECTOR_DEGREE = 2; /*! @todo: Expose to ParameterHandler */

const double LIQUID_DYNAMIC_VISOCITY = 1.; /*! @todo: Expose to ParameterHandler */

const double EPSILON = 1.e-14; /*! @todo: Expose to ParameterHandler */

const bool WRITE_LINEAR_SYSTEM = true; /*! @todo: Expose to ParameterHandler */

#endif