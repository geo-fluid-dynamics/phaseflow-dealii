#ifndef my_grid_generator_h
#define my_grid_generator_h

#include <iostream>
#include <cmath>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

namespace MyGridGenerator
{
    using namespace dealii;
    
    template<int dim>
    void create_coarse_grid(
        Triangulation<dim> &triangulation,
        unsigned int &boundary_count,
        const std::string grid_name,
        const std::vector<double> sizes)
    {
        if (grid_name == "hyper_cube")
        {
            GridGenerator::hyper_cube(
                triangulation,
                sizes[0],
                sizes[1],
                true);

            boundary_count = pow(2, dim);
        }
        else if (grid_name == "hyper_rectangle")
        {
            GridGenerator::hyper_rectangle(
                triangulation,
                {sizes[0], sizes[1]},
                {sizes[2], sizes[3]},
                true);

            boundary_count = pow(2, dim);
        }
        else
        {
            throw(ExcNotImplemented());
        }
    }
 
}

#endif
