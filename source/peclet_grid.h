#include "refinement.h"

template<>
void Peclet<1>::create_coarse_grid()
{
    const unsigned int dim = 1;
    // Create grid
    MyGridGenerator::create_coarse_grid(
        this->triangulation,
        this->manifold_ids, this->manifold_descriptors,
        this->params.geometry.grid_name, params.geometry.sizes);
    // Shift and rotate the grid.
    Point<dim> shifted_center;
    for (unsigned int i = 0; i < dim; i++)
    {
        shifted_center[i] = params.geometry.transformations[i];
    }
    GridTools::shift(shifted_center, this->triangulation); 
    this->spherical_manifold_center = shifted_center;
}

template <>
void Peclet<2>::create_coarse_grid()
{
    const unsigned int dim = 2;
    // Create grid
    MyGridGenerator::create_coarse_grid(
        this->triangulation,
        this->manifold_ids, this->manifold_descriptors,
        params.geometry.grid_name, params.geometry.sizes);
    // Shift and rotate the grid.
    Point<dim> shifted_center;
    GridTools::rotate(params.geometry.transformations[2], this->triangulation);
    for (unsigned int i = 0; i < dim; i++)
    {
        shifted_center[i] = params.geometry.transformations[i];
    }
    GridTools::shift(shifted_center, this->triangulation); 
    this->spherical_manifold_center = shifted_center;
}

template<>
void Peclet<3>::create_coarse_grid()
{
    Assert(false, ExcNotImplemented()); // Only missing a 3D rotation method
}

template<int dim>
void Peclet<dim>::adaptive_refine()
{
    SolutionTransfer<dim> solution_trans(this->dof_handler);
    Vector<double> previous_solution;
    previous_solution = this->solution;
    Refinement::adaptive_refine_mesh<dim>();
    this->setup_system();
    solution_trans.interpolate(previous_solution, this->solution);
    this->constraints.distribute(this->solution);
}
  