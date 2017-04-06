#include "refinement.h"

template <int dim>
void Peclet<dim>::create_coarse_grid()
{
    MyGridGenerator::create_coarse_grid(
        this->triangulation,
        this->boundary_count,
        this->params.geometry.grid_name, params.geometry.sizes);
}
