#ifndef peclet_global_parameters_h
#define peclet_global_parameters_h

/*! Set global parameters.

Most of these aren't actually required at compile time, and should ideally be exposed to ParameterHandler.
This is just a temporary solution.

*/
const bool ENERGY_ENABLED = true; /*! @todo: Expose to ParameterHandler */

const double RAYLEIGH_NUMBER = 1.e6; /*! @todo: Expose to ParameterHandler */

const double PRANDTL_NUMBER = 0.71; /*! @todo: Expose to ParameterHandler */

const double REYNOLDS_NUMBER = 1.;

const double LIQUID_DYNAMIC_VISOCITY = 1.; /*! @todo: Expose to ParameterHandler */

const double SOLID_CONDUCTIVITY = 1.; /*! @todo: Expose to ParameterHandler */

const double LIQUID_CONDUCTIVITY = 2.; /*! @todo: Expose to ParameterHandler */

const unsigned int SCALAR_DEGREE = 1; /*! @todo: Expose to ParameterHandler */

const double EPSILON = 1.e-14; /*! @todo: Expose to ParameterHandler */

const bool WRITE_LINEAR_SYSTEM = true; /*! @todo: Expose to ParameterHandler */

const double STEADY_TOLERANCE = 1.e-4; /*! @todo: Expose to ParameterHandler */

const double TIME_GROWTH_RATE = 2.; /*! @todo: Expose to ParameterHandler */

const std::set<dealii::types::boundary_id> ADIABATIC_WALLS = {2, 3}; /*! @todo: Generalize boundary conditions */

#endif
