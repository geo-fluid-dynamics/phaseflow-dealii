#include "refinement.h"

template <int dim>
void Peclet<dim>::create_coarse_grid()
{
    MyGridGenerator::create_coarse_grid(
        this->triangulation,
        this->manifold_ids,
        this->manifold_descriptors,
        this->boundary_count,
        this->params.geometry.grid_name,
        params.geometry.sizes);
}
